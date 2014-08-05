/*
 * Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2006, Technische Universitaet Muenchen, and individual
 * contributors as indicated by the @authors tag. See the 
 * copyright.txt in the distribution for a full listing of individual
 * contributors.
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 */

/**
 * @ingroup tracking_algorithms
 * @file
 * Implements functions for 2D-3D pose estimation.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

 
#include "PlanarPoseEstimation.h"
#include "../Function/MultiplePointProjection.h"
#include "../Function/MultiplePointProjectionError.h"
#include "../Function/MultipleCameraProjectionError.h"
#include "../Function/ProjectivePoseNormalize.h"

#include "../Homography.h"
#include "../Projection.h"


#include <utMath/VectorFunctions.h>
#include <utMath/MatrixOperations.h>
#include <utMath/Stochastic/BackwardPropagation.h>


#include <math.h>
#include <iostream>


#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>


#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#endif

//#define OPTIMIZATION_LOGGING
// #include <log4cpp/Category.hh>
// static log4cpp::Category& optLogger( log4cpp::Category::getInstance( "Ubitrack.Algorithm.2D3DPoseEstimation" ) );
#include <utMath/Optimization/LevenbergMarquardt.h>


// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;
using namespace Ubitrack::Math;

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation2D3D {

/** \internal */
template< typename T >
Math::Pose poseFromHomographyImpl( const Matrix< T, 3, 3 >& H, const Matrix< T, 3, 3 >& invK )
{
	// compute R = K^-1 H
	Matrix< T, 3, 3 > R( ublas::prod( invK, H ) );
	
	// CW@2013-12-03 changed this lately from negative to positive
	// due to serious problems in the estimation of poses.
	// still do not know what is correct: positive or negative?
	// all my tests seem to point to positive, therefore :
	// make sure the z-coordinate is positive
	if ( R( 2, 2 ) > 0 )
		R *= -1.0;

	// compute length of the first two colums
	const T fXLen = ublas::norm_2( ublas::column( R, 0 ) );
	const T fYLen = ublas::norm_2( ublas::column( R, 1 ) );
	
	// copy & normalize translation
	const T fTransScale = 2.0f / ( fXLen + fYLen );
	Vector< T, 3 > t( ublas::column( R, 2 ) * fTransScale );
		
#if defined( HAVE_LAPACK ) && !defined( __APPLE__ ) // APPLE vecLib sucks

	// perform svd-based orthogonalization
	Math::Matrix< T, 3, 3 > u;
	Math::Matrix< T, 3, 3 > right;
	Math::Vector< T, 2 > s;
	//CW@2013-12-06:
	// last change was wrong, s needs to be set to 2
	// actually the svd looks like: R_3x2 * S_2 * Vt_2x2, although matrices are 3x3
	ublas::matrix_range< Math::Matrix< T, 3, 3 > > Rleft( R, ublas::range( 0, 3 ), ublas::range( 0, 2 ) );
	ublas::matrix_range< Math::Matrix< T, 3, 3 > > vt( right, ublas::range( 0, 2 ), ublas::range( 0, 2 ) );
	boost::numeric::bindings::lapack::gesvd( 'A', 'A', Rleft, s, u, vt );
	
	right( 0, 2 ) = right( 1, 2 ) = 0;
	right( 2, 0 ) = right( 2, 1 ) = 0;
	right( 2, 2 ) = Math::determinant( vt ) * Math::determinant( u ); // should be -1 or +1
	R = ublas::prod( u, right );

#else

	// normalize first two colums
	ublas::column( R, 0 ) /= fXLen;
	ublas::column( R, 1 ) /= fYLen;

	// compute third row as vector product
	ublas::column( R, 2 ) = cross_product( ublas::column( R, 0 ), ublas::column( R, 1 ) );
	
	// normalize cross product
	const T fZLen = ublas::norm_2( ublas::column( R, 2 ) );
	ublas::column( R, 2 ) /= fZLen;
	
	// recompute y vector from x and z
	ublas::column( R, 1 ) = cross_product( ublas::column( R, 2 ), ublas::column( R, 0 ) );

#endif
	// compute rotation quaternion from matrix
	return Pose( Quaternion( R ), t );
}

Math::Pose poseFromHomography( const Math::Matrix< float, 3, 3 >& H, const Math::Matrix< float, 3, 3 >& invK )
{
	return poseFromHomographyImpl( H, invK );
}

Math::Pose poseFromHomography( const Math::Matrix< double, 3, 3 >& H, const Math::Matrix< double, 3, 3 >& invK )
{
	return poseFromHomographyImpl( H, invK );
}


#ifdef HAVE_LAPACK

/** \internal */	
template< typename T > 
T optimizePoseImpl( Pose& p, const std::vector< Vector< T, 2 > >& p2D, const std::vector< Vector< T, 3 > >& p3D, 
	const Matrix< T, 3, 3 >& cam, const std::size_t nIterations  )
{
	// copy rot & trans to parameter vector
	Vector< T, 7 > params;
	p.toVector( params );

	// copy 2D points to measurement vector
	Math::Vector< T >  measurements( 2 * p2D.size() );
	for ( std::size_t i( 0 ); i < p2D.size(); i++ )
		ublas::subrange( measurements, 2*i, (i+1)*2 ) = p2D[ i ];

	// perform optimization
	Function::MultiplePointProjection< T > projection( p3D, cam );
	T fRes = Optimization::levenbergMarquardt( projection, params, measurements, 
		Optimization::OptTerminate( nIterations, 1e-6 ), Function::ProjectivePoseNormalize() );

	// copy back rot & trans from vector
	p = Pose::fromVector( params );
		
	return fRes;
}

float optimizePose( Math::Pose& p, const std::vector< Math::Vector< float, 2 > >& p2D, 
	const std::vector< Math::Vector< float, 3 > >& p3D, const Math::Matrix< float, 3, 3 >& cam,
	const std::size_t nIterations )
{
	return optimizePoseImpl( p, p2D, p3D, cam, nIterations );
}

double optimizePose( Math::Pose& p, const std::vector< Math::Vector< double, 2 > >& p2D, 
	const std::vector< Math::Vector< double, 3 > >& p3D, const Math::Matrix< double, 3, 3 >& cam,
	const std::size_t nIterations )
{
	return optimizePoseImpl( p, p2D, p3D, cam, nIterations );
}


/** \internal */
template< typename T >
Matrix< T, 6, 6 > singleCameraPoseErrorImpl( const Pose& p, const std::vector< Vector< T, 3 > >& p3D, 
	const Matrix< T, 3, 3 >& cam, T imageError )
{
	// copy rot & trans to parameter vector
	Vector< T, 7 > params;
	p.toVector( params );

	// calculate error
	Matrix< T, 6, 6 > result;
	Function::MultiplePointProjectionError< T > projection( p3D, cam );

	Stochastic::backwardPropagationIdentity( result, imageError, projection, params );
	
	return result;
}

Matrix< float, 6, 6 > singleCameraPoseError( const Math::Pose& p, const std::vector< Math::Vector< float, 3 > >& p3D, 
	const Math::Matrix< float, 3, 3 >& cam, float imageError )
{
	return singleCameraPoseErrorImpl( p, p3D, cam, imageError );
}

Matrix< double, 6, 6 > singleCameraPoseError( const Math::Pose& p, const std::vector< Math::Vector< double, 3 > >& p3D, 
	const Math::Matrix< double, 3, 3 >& cam, double imageError )
{
	return singleCameraPoseErrorImpl( p, p3D, cam, imageError );
}


/** \internal */
template< typename T >
Math::Matrix< T, 6, 6 > multipleCameraPoseErrorImpl( const Math::Pose& p, 
	const std::vector< Math::Vector< T, 3 > >& p3D, 
	const std::vector< Math::Matrix< T, 3, 4 > >& cameras, 
	const std::vector< std::pair< std::size_t, std::size_t > > observations, 
	T imageError )
{
	// copy rot & trans to parameter vector
	Vector< T, 7 > params;
	p.toVector( params );

	// calculate error
	Matrix< T, 6, 6 > result;
	Function::MultipleCameraProjectionError< T > projection( p3D, cameras, observations );

	Stochastic::backwardPropagationIdentity( result, imageError, projection, params );
	
	return result;
}
	
Math::Matrix< float, 6, 6 > multipleCameraPoseError( const Math::Pose& p, 
	const std::vector< Math::Vector< float, 3 > >& p3D, 
	const std::vector< Math::Matrix< float, 3, 4 > >& cameras, 
	const std::vector< std::pair< std::size_t, std::size_t > > observations, 
	float imageError )
{
	return multipleCameraPoseErrorImpl( p, p3D, cameras, observations, imageError );
}
	
Math::Matrix< double, 6, 6 > multipleCameraPoseError( const Math::Pose& p, 
	const std::vector< Math::Vector< double, 3 > >& p3D, 
	const std::vector< Math::Matrix< double, 3, 4 > >& cameras, 
	const std::vector< std::pair< std::size_t, std::size_t > > observations, 
	double imageError )
{
	return multipleCameraPoseErrorImpl( p, p3D, cameras, observations, imageError );
}

double reprojectionError( 
	const std::vector< Math::Vector< double, 2 > >& p2d,
	const std::vector< Math::Vector< double, 3 > >& p3d, 
	Math::Pose p, 
	Math::Matrix< double, 3, 3 > cam )
{
	Math::Matrix< double, 3, 4 > projMat;
	double res = 0.0;

	// Create a pose matrix
	Math::Matrix< double, 3, 3 > rot = p.rotation();
	Math::Vector< double, 3 > trans = p.translation();
	projMat(0,0) = rot(0,0);
	projMat(0,1) = rot(0,1);
	projMat(0,2) = rot(0,2);
	projMat(1,0) = rot(1,0);
	projMat(1,1) = rot(1,1);
	projMat(1,2) = rot(1,2);
	projMat(2,0) = rot(2,0);
	projMat(2,1) = rot(2,1);
	projMat(2,2) = rot(2,2);
	projMat(0,3) = trans[0];
	projMat(1,3) = trans[1];
	projMat(2,3) = trans[2];
	
	// P = K * [R|T] 
	projMat = boost::numeric::ublas::prod(cam, projMat) ; 
	
	// 2D = P * 3D Point
	// Reproject 3D-points to 2D points.
	const std::size_t n_points( p3d.size() );
	for( std::size_t i( 0 ); i < n_points; i++) 
	{	
		Math::Vector< double, 4 > hom( (p3d.at(i))[0], (p3d.at(i))[1], (p3d.at(i))[2], 1.0 );
		Math::Vector< double, 3 > tmp;
		
		tmp = boost::numeric::ublas::prod(projMat, hom) ; 
		double w = tmp[2];
		tmp = tmp / w;
		
		res += ( tmp[0] - (p2d.at(i))[0] ) * ( tmp[0] - (p2d.at(i))[0] ) + ( tmp[1] - (p2d.at(i))[1] ) * ( tmp[1] - (p2d.at(i))[1] );
	}
	/*
	res /= p3d.size();
	res = sqrt( res );
	*/
	return res;
}

Math::ErrorPose computePose( 
		const std::vector< Math::Vector< double, 2 > >& p2d,
		const std::vector< Math::Vector< double, 3 > >& p3d,
		const Math::Matrix< double, 3, 3 >& cam,
		bool optimize,
		enum InitializationMethod initMethod
	)
{
	double dummy;
	return computePose( p2d, p3d, cam, dummy, optimize, initMethod );
}

Math::ErrorPose computePose( 
		const std::vector< Math::Vector< double, 2 > >& p2d,
		const std::vector< Math::Vector< double, 3 > >& p3d,
		const Math::Matrix< double, 3, 3 >& cam,
		double& residual,
		bool optimize,
		enum InitializationMethod initMethod
	)
{
	const std::size_t n_points( p2d.size() );
	if ( n_points < 4 ) {
		UBITRACK_THROW( "2D3D pose estimation configured to use at least 4 points" );
	}
	
	OPT_LOG_DEBUG( "Performing pose estimation using " << n_points << " points" );
	OPT_LOG_TRACE( "2D points: " << p2d );
	OPT_LOG_TRACE( "3D points: " << p3d );

	// invert camera matrix
	Math::Matrix< double, 3, 3 > invK( Math::invert_matrix( cam ) );

	Math::Pose pose;

	bool bInitialized = false;
	
	if ( initMethod == NONPLANAR_PROJECTION && n_points >= 6 )
	{
		// initialize from 3x4 projection matrix
		Math::Matrix< double, 3, 4 > P( Algorithm::projectionDLT( p3d, p2d ) );
		Math::Matrix< double, 3, 4 > Rt( ublas::prod( invK, P ) );
		OPT_LOG_TRACE( "inital [R|t]: " << std::endl << Rt );

		//TODO ### Just a check
		Math::Matrix< double, 3, 3 > kTest;
		Math::Matrix< double, 3, 3 > rTest;
		Math::Vector< double, 3 > tTest;
		decomposeProjection( kTest, rTest, tTest, P ); 
		OPT_LOG_TRACE( "K (given): " << std::endl << cam );
		OPT_LOG_TRACE( "K from decomposition of P: " << std::endl << kTest );
		OPT_LOG_TRACE( "R from decomposition of P: " << std::endl << rTest );
		OPT_LOG_TRACE( "t from decomposition of P: " << std::endl << tTest );

		// perform svd decomposition to get a pure rotation matrix
		Math::Matrix< double, 3, 3 > u;
		Math::Matrix< double, 3, 3 > vt;
		Math::Vector< double, 3 > s;
		ublas::matrix_range< Math::Matrix< double, 3, 4 > > R( Rt, ublas::range( 0, 3 ), ublas::range( 0, 3 ) );
		ublas::matrix_column< Math::Matrix< double, 3, 4 > > t( Rt, 3 );

		if ( Math::determinant( R ) < 0 )
			Rt *= -1;

		boost::numeric::bindings::lapack::gesvd( 'A', 'A', R, s, u, vt );
		OPT_LOG_TRACE( "s: " << s );
		OPT_LOG_TRACE( "U: " << std::endl << u );
		OPT_LOG_TRACE( "V^T: " << std::endl << vt );

		// Compute condition number to check the "orthonormality" of the obtained rotation matrix
		if ( ( s( 0 ) / s( 2 ) ) < 2 )
		{
			R = ublas::prod( u, vt );
			t /= s( 0 ) * s( 1 ) * s( 2 ); // det( original-R )

			pose = Math::Pose( Math::Quaternion( R ), t );
			bInitialized = true;

			OPT_LOG_TRACE( "Pose from projection matrix: " << pose );
		}
		else 
		{
			OPT_LOG_DEBUG( "3x4 DLT was unstable (planar target?)" );
		}
	}

	if ( !bInitialized ) 
	{
	
		// 1st possibility:
		// first points start with same z-value (coplanar). Take all 
		// points in a row that have the same dimension and calculate
		// a first homography of these points if there are more than four.
		// this assumption works in cases of markers (4 points) and planar 
		// calibration grid structures with much more than 4 points.
		std::vector< Math::Vector< double, 2 > > p3dAs2d;
		p3dAs2d.reserve( n_points );
		double last_dim = p3d[ 0 ]( 2 );
			
		for ( std::size_t i( 0 ); i < n_points; ++i )
		{
			if( p3d[ i ]( 2 ) == last_dim )
				p3dAs2d.push_back( Math::Vector< double, 2 >( p3d[ i ]( 0 ), p3d[ i ]( 1 ) ) );
			else
				break;
		}
		const std::size_t n3D = p3dAs2d.size();
			
		if( n3D > 3 )
		{
			// copy corresponding elements to a temporary vector
			// and compute the homography
			Math::Matrix< double, 3, 3 > H =  Algorithm::homographyDLT( p3dAs2d, std::vector< Math::Vector< double, 2 > >( p2d.begin(), p2d.begin() + n3D ) );
			
			OPT_LOG_TRACE( "Homography: " << H );

			// compute initial pose from homography
			pose = Algorithm::PoseEstimation2D3D::poseFromHomography( H, invK );
			
			OPT_LOG_TRACE( "Pose from homography: " << pose )
		}
		else
		// 2nd possibility:
		// first points lye within a rotated plane which is not parallel
		// to the xy-plane. Use the first four values to calculate 
		// an initial homography and a corresponding pose.
		{
			// compute a rotation matrix that will bring the points into a plane with equal z
			Math::Vector< double, 3 > vX( p3d[ 1 ] - p3d[ 0 ] );
			double f = ublas::norm_2( vX );
			vX /= f;

			Math::Vector< double, 3 > vZ( p3d[ 2 ] - p3d[ 0 ] );
			f = ublas::norm_2( vZ );
			vZ /= f;

			// Check whether first three points are colinear
			double collinearity = inner_prod( vX, vZ ) > 0.8;
			OPT_LOG_TRACE( "Checking collinearity constraint (should be lower than 0.8): " << collinearity );
			if ( collinearity ) {
				OPT_LOG_TRACE( "Points are collinear" );
				UBITRACK_THROW( "Pose estimation requires four coplanar points in general position but three of them are collinear" );
			}

			vZ = Math::cross_product( vX, vZ );
			f = ublas::norm_2( vZ );
			vZ /= f;

			Math::Matrix< double, 3, 3 > P;
			ublas::row( P, 0 ) = vX;
			ublas::row( P, 2 ) = vZ;
			ublas::row( P, 1 ) = Math::cross_product( vZ, vX );

			// compute a translation
			Math::Vector< double, 3 > t( -ublas::prod( P, p3d[ 0 ] ) );

			OPT_LOG_TRACE( "Computed alignment, now checking coplanarity constraint..." );

			p3dAs2d.clear();
			for ( std::size_t i( 0 ); i < 4; i++ )
			{
				Math::Vector< double, 3 > p3dtrans = ublas::prod( P, p3d[ i ] ) + t;
				OPT_LOG_TRACE( "z-value of point " << i << ": " << fabs( p3dtrans( 2 ) ) );
				if ( fabs( p3dtrans( 2 ) ) > 1e-2 ) {
					OPT_LOG_TRACE( "Points are NOT very coplanar" );
					//TODO ### UBITRACK_THROW( "Pose estimation requires four coplanar points" );
				}
				p3dAs2d.push_back( Math::Vector< double, 2 >( p3dtrans( 0 ), p3dtrans( 1 ) ) );
			}

			// compute homography
			Math::Matrix< double, 3, 3 > H;
			if ( n_points > 4 )
				// copy first four elements to new vector
				H =  Algorithm::homographyDLT( p3dAs2d, std::vector< Math::Vector< double, 2 > >( p2d.begin(), p2d.begin() + 4 ) );
			else
				H =  Algorithm::homographyDLT( p3dAs2d, p2d );

			// compute initial pose from homography
			pose = Algorithm::PoseEstimation2D3D::poseFromHomography( H, invK ) * Math::Pose( Math::Quaternion( P ), t );
			
			OPT_LOG_TRACE( "Pose from homography (rotated): " << pose );
			Math::Matrix< double, 3, 3 > rotMat;
			pose.rotation().toMatrix( rotMat );
			OPT_LOG_TRACE( "Rotation matrix (rotated): " << rotMat );
		}
	}
	
	// non-linear minimization
	Math::Matrix< double, 6, 6 > covMatrix;
	if ( optimize )
	{
		residual = Algorithm::PoseEstimation2D3D::optimizePose( pose, p2d, p3d, cam );
		OPT_LOG_DEBUG( "Refined pose: " << pose << ", residual of 2D image measurements: " << residual);
	}
	else
	{
		residual = reprojectionError( p2d, p3d, pose, cam );
		OPT_LOG_DEBUG( "NOT refined pose: " << pose << ", residual of 2D image measurements: " << residual);	
	}
	
	covMatrix = Algorithm::PoseEstimation2D3D::singleCameraPoseError( pose, p3d, cam, residual );	
	residual = sqrt( residual / ( n_points * 2 ) );
	
	return Math::ErrorPose( pose, covMatrix );
}

#endif // HAVE_LAPACK

} } } // namespace Ubitrack::Algorithm::PoseEstimation2D3D

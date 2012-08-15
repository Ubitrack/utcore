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
 * Implementats functions for 2D-3D pose estimation.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#include <iostream>
#include <math.h>

#include <log4cpp/Category.hh>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <utMath/MatrixOperations.h>
#include <utMath/BackwardPropagation.h>
#include <utCalibration/Homography.h>
#include <utCalibration/Projection.h>


#include "2D3DPoseEstimation.h"
#include "Function/MultiplePointProjection.h"
#include "Function/MultiplePointProjectionError.h"
#include "Function/MultipleCameraProjectionError.h"
#include "Function/ProjectivePoseNormalize.h"

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#endif

#define OPTIMIZATION_LOGGING
static log4cpp::Category& optLogger( log4cpp::Category::getInstance( "Ubitrack.Calibration.2D3DPoseEstimation" ) );
#include <utMath/LevenbergMarquardt.h>


// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;
using namespace Ubitrack::Math;

namespace Ubitrack { namespace Calibration {

/** \internal */
template< typename T >
Math::Pose poseFromHomographyImpl( const Matrix< 3, 3, T > H, const Matrix< 3, 3, T > invK )
{
	// compute R = K^-1 H
	Matrix< 3, 3 > R( ublas::prod( invK, H ) );
	
	// make sure the z-coordinate is negative
	if ( R( 2, 2 ) > 0 )
		R *= -1.0;

	// compute length of the first two colums
	double fXLen = ublas::norm_2( ublas::column( R, 0 ) );
	double fYLen = ublas::norm_2( ublas::column( R, 1 ) );
	
	// copy & normalize translation
	double fTransScale = 2.0f / ( fXLen + fYLen );
	Vector< 3 > t( ublas::column( R, 2 ) * fTransScale );
		
#if defined( HAVE_LAPACK ) && !defined( __APPLE__ ) // APPLE vecLib sucks

	// perform svd-based orthogonalization
	Math::Matrix< 3, 3 > u;
	Math::Matrix< 3, 3 > right;
	Math::Vector< 2 > s;
	ublas::matrix_range< Matrix< 3, 3 > > Rleft( R, ublas::range( 0, 3 ), ublas::range( 0, 2 ) );
	ublas::matrix_range< Matrix< 3, 3 > > vt( right, ublas::range( 0, 2 ), ublas::range( 0, 2 ) );
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
	ublas::column( R, 2 ) = cross_prod( ublas::column( R, 0 ), ublas::column( R, 1 ) );
	
	// normalize cross product
	double fZLen = ublas::norm_2( ublas::column( R, 2 ) );
	ublas::column( R, 2 ) /= fZLen;
	
	// recompute y vector from x and z
	ublas::column( R, 1 ) = cross_prod( ublas::column( R, 2 ), ublas::column( R, 0 ) );

#endif
	// compute rotation quaternion from matrix
	return Pose( Quaternion( R ), t );
}

Math::Pose poseFromHomography( const Math::Matrix< 3, 3, float > H, const Math::Matrix< 3, 3, float > invK )
{
	return poseFromHomographyImpl( H, invK );
}

Math::Pose poseFromHomography( const Math::Matrix< 3, 3, double > H, const Math::Matrix< 3, 3, double > invK )
{
	return poseFromHomographyImpl( H, invK );
}


#ifdef HAVE_LAPACK

/** \internal */	
template< typename T > 
T optimizePoseImpl( Pose& p, const std::vector< Vector< 2, T > >& p2D, const std::vector< Vector< 3, T > >& p3D, 
	const Matrix< 3, 3, T >& cam, unsigned nIterations  )
{
	// copy rot & trans to parameter vector
	Vector< 7, T > params;
	p.toVector( params );

	// copy 2D points to measurement vector
	ublas::vector< T >  measurements( 2 * p2D.size() );
	for ( unsigned i = 0; i < p2D.size(); i++ )
		ublas::subrange( measurements, 2*i, (i+1)*2 ) = p2D[ i ];

	// perform optimization
	Function::MultiplePointProjection< T > projection( p3D, cam );
	T fRes = levenbergMarquardt( projection, params, measurements, 
		OptTerminate( nIterations, 1e-6 ), Function::ProjectivePoseNormalize() );

	// copy back rot & trans from vector
	p = Pose::fromVector( params );
		
	return fRes;
}

float optimizePose( Math::Pose& p, const std::vector< Math::Vector< 2, float > >& p2D, 
	const std::vector< Math::Vector< 3, float > >& p3D, const Math::Matrix< 3, 3, float >& cam,
	unsigned nIterations )
{
	return optimizePoseImpl( p, p2D, p3D, cam, nIterations );
}

double optimizePose( Math::Pose& p, const std::vector< Math::Vector< 2, double > >& p2D, 
	const std::vector< Math::Vector< 3, double > >& p3D, const Math::Matrix< 3, 3, double >& cam,
	unsigned nIterations )
{
	return optimizePoseImpl( p, p2D, p3D, cam, nIterations );
}


/** \internal */
template< typename T >
Matrix< 6, 6, T > singleCameraPoseErrorImpl( const Pose& p, const std::vector< Vector< 3, T > >& p3D, 
	const Matrix< 3, 3, T >& cam, T imageError )
{
	// copy rot & trans to parameter vector
	Vector< 7, T > params;
	p.toVector( params );

	// calculate error
	Matrix< 6, 6, T > result;
	Function::MultiplePointProjectionError< T > projection( p3D, cam );

	backwardPropagationIdentity( result, imageError, projection, params );
	
	return result;
}

Matrix< 6, 6, float > singleCameraPoseError( const Math::Pose& p, const std::vector< Math::Vector< 3, float > >& p3D, 
	const Math::Matrix< 3, 3, float >& cam, float imageError )
{
	return singleCameraPoseErrorImpl( p, p3D, cam, imageError );
}

Matrix< 6, 6, double > singleCameraPoseError( const Math::Pose& p, const std::vector< Math::Vector< 3, double > >& p3D, 
	const Math::Matrix< 3, 3, double >& cam, double imageError )
{
	return singleCameraPoseErrorImpl( p, p3D, cam, imageError );
}


/** \internal */
template< typename T >
Math::Matrix< 6, 6, T > multipleCameraPoseErrorImpl( const Math::Pose& p, 
	const std::vector< Math::Vector< 3, T > >& p3D, 
	const std::vector< Math::Matrix< 3, 4, T > >& cameras, 
	const std::vector< std::pair< unsigned, unsigned > > observations, 
	T imageError )
{
	// copy rot & trans to parameter vector
	Vector< 7, T > params;
	p.toVector( params );

	// calculate error
	Matrix< 6, 6, T > result;
	Function::MultipleCameraProjectionError< T > projection( p3D, cameras, observations );

	backwardPropagationIdentity( result, imageError, projection, params );
	
	return result;
}
	
Math::Matrix< 6, 6, float > multipleCameraPoseError( const Math::Pose& p, 
	const std::vector< Math::Vector< 3, float > >& p3D, 
	const std::vector< Math::Matrix< 3, 4, float > >& cameras, 
	const std::vector< std::pair< unsigned, unsigned > > observations, 
	float imageError )
{
	return multipleCameraPoseErrorImpl( p, p3D, cameras, observations, imageError );
}
	
Math::Matrix< 6, 6, double > multipleCameraPoseError( const Math::Pose& p, 
	const std::vector< Math::Vector< 3, double > >& p3D, 
	const std::vector< Math::Matrix< 3, 4, double > >& cameras, 
	const std::vector< std::pair< unsigned, unsigned > > observations, 
	double imageError )
{
	return multipleCameraPoseErrorImpl( p, p3D, cameras, observations, imageError );
}

double reprojectionError( 
	const std::vector< Math::Vector< 2 > >& p2d,
	const std::vector< Math::Vector< 3 > >& p3d, 
	Math::Pose p, 
	Math::Matrix< 3, 3 > cam )
{
	Math::Matrix< 3, 4 > projMat;
	double res = 0.0;

	// Create a pose matrix
	Math::Matrix< 3, 3 > rot = p.rotation();
	Math::Vector< 3 > trans = p.translation();
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
	for( unsigned int i = 0; i < p3d.size(); i++) 
	{	
		Math::Vector< 4 , double> hom( (p3d.at(i))[0], (p3d.at(i))[1], (p3d.at(i))[2], 1.0 );
		Math::Vector< 3 , double> tmp;
		
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
		const std::vector< Math::Vector< 2 > >& p2d,
		const std::vector< Math::Vector< 3 > >& p3d,
		const Math::Matrix< 3, 3 >& cam,
		bool optimize,
		enum InitializationMethod initMethod
	)
{
	double dummy;
	return computePose( p2d, p3d, cam, dummy, optimize, initMethod );
}

Math::ErrorPose computePose( 
		const std::vector< Math::Vector< 2 > >& p2d,
		const std::vector< Math::Vector< 3 > >& p3d,
		const Math::Matrix< 3, 3 >& cam,
		double& residual,
		bool optimize,
		enum InitializationMethod initMethod
	)
{
	if ( p2d.size() < 4 ) {
		UBITRACK_THROW( "2D3D pose estimation configured to use at least 4 points" );
	}
	
	LOG4CPP_DEBUG( optLogger, "Performing pose estimation using " << p2d.size() << " points" );
	LOG4CPP_TRACE( optLogger, "2D points: " << p2d );
	LOG4CPP_TRACE( optLogger, "3D points: " << p3d );

	// invert camera matrix
	Math::Matrix< 3, 3 > invK( Math::invert_matrix( cam ) );

	Math::Pose pose;

	bool bInitialized = false;
	
	if ( initMethod == NONPLANAR_PROJECTION && p2d.size() >= 6 )
	{
		// initialize from 3x4 projection matrix
		Math::Matrix< 3, 4 > P( Calibration::projectionDLT( p3d, p2d ) );
		Math::Matrix< 3, 4 > Rt( ublas::prod( invK, P ) );
		LOG4CPP_TRACE( optLogger, "inital [R|t]: " << std::endl << Rt );

		//TODO ### Just a check
		Math::Matrix< 3, 3, float > kTest;
		Math::Matrix< 3, 3, float > rTest;
		Math::Vector< 3, float > tTest;
		decomposeProjection( kTest, rTest, tTest, P ); 
		LOG4CPP_TRACE( optLogger, "K (given): " << std::endl << cam );
		LOG4CPP_TRACE( optLogger, "K from decomposition of P: " << std::endl << kTest );
		LOG4CPP_TRACE( optLogger, "R from decomposition of P: " << std::endl << rTest );
		LOG4CPP_TRACE( optLogger, "t from decomposition of P: " << std::endl << tTest );

		// perform svd decomposition to get a pure rotation matrix
		Math::Matrix< 3, 3 > u;
		Math::Matrix< 3, 3 > vt;
		Math::Vector< 3 > s;
		ublas::matrix_range< Math::Matrix< 3, 4 > > R( Rt, ublas::range( 0, 3 ), ublas::range( 0, 3 ) );
		ublas::matrix_column< Math::Matrix< 3, 4 > > t( Rt, 3 );

		if ( Math::determinant( R ) < 0 )
			Rt *= -1;

		boost::numeric::bindings::lapack::gesvd( 'A', 'A', R, s, u, vt );
		LOG4CPP_TRACE( optLogger, "s: " << s );
		LOG4CPP_TRACE( optLogger, "U: " << std::endl << u );
		LOG4CPP_TRACE( optLogger, "V^T: " << std::endl << vt );

		// Compute condition number to check the "orthonormality" of the obtained rotation matrix
		if ( ( s( 0 ) / s( 2 ) ) < 2 )
		{
			R = ublas::prod( u, vt );
			t /= s( 0 ) * s( 1 ) * s( 2 ); // det( original-R )

			pose = Math::Pose( Math::Quaternion( R ), t );
			bInitialized = true;

			LOG4CPP_TRACE( optLogger, "Pose from projection matrix: " << pose );
		}
		else 
		{
			LOG4CPP_DEBUG( optLogger, "3x4 DLT was unstable (planar target?)" );
		}
	}

	if ( !bInitialized ) 
	{
		if ( p3d[ 0 ]( 2 ) == 0 && p3d[ 1 ]( 2 ) == 0 && p3d[ 2 ]( 2 ) == 0 && p3d[ 3 ]( 2 ) == 0 )
		{
			// markers already have z=0
			std::vector< Math::Vector< 2 > > p3dAs2d;
			for ( unsigned i = 0; i < 4; i++ )
				p3dAs2d.push_back( Math::Vector< 2 >( p3d[ i ]( 0 ), p3d[ i ]( 1 ) ) );

			// compute homography
			Math::Matrix< 3, 3 > H;
			if ( p2d.size() > 4 )
				// copy first four elements to new vector
				H =  Calibration::homographyDLT( p3dAs2d, std::vector< Math::Vector< 2 > >( p2d.begin(), p2d.begin() + 4 ) );
			else
				H =  Calibration::homographyDLT( p3dAs2d, p2d );
			LOG4CPP_TRACE( optLogger, "Homography: " << H );

			// compute initial pose from homography
			pose = Calibration::poseFromHomography( H, invK );
			LOG4CPP_TRACE( optLogger, "Pose from homography: " << pose );
		}
		else
		{
			// compute a rotation matrix that will bring the points into a plane with equal z
			Math::Vector< 3 > vX( p3d[ 1 ] - p3d[ 0 ] );
			double f = ublas::norm_2( vX );
			vX /= f;

			Math::Vector< 3 > vZ( p3d[ 2 ] - p3d[ 0 ] );
			f = ublas::norm_2( vZ );
			vZ /= f;

			// Check whether first three points are colinear
			double colinearity = inner_prod( vX, vZ ) > 0.8;
			LOG4CPP_TRACE( optLogger, "Checking colinearity constraint (should be lower than 0.8): " << colinearity );
			if ( colinearity ) {
				LOG4CPP_TRACE( optLogger, "Points are colinear" );
				UBITRACK_THROW( "Pose estimation requires four coplanar points in general position but three of them are colinear" );
			}

			vZ = Math::cross_prod( vX, vZ );
			f = ublas::norm_2( vZ );
			vZ /= f;

			Math::Matrix< 3, 3 > P;
			ublas::row( P, 0 ) = vX;
			ublas::row( P, 2 ) = vZ;
			ublas::row( P, 1 ) = Math::cross_prod( vZ, vX );

			// compute a translation
			Math::Vector< 3 > t( -ublas::prod( P, p3d[ 0 ] ) );

			LOG4CPP_TRACE( optLogger, "Computed alignment, now checking coplanarity constraint..." );

			std::vector< Math::Vector< 2 > > p3dAs2d;
			for ( unsigned i = 0; i < 4; i++ )
			{
				Math::Vector< 3 > p3dtrans = ublas::prod( P, p3d[ i ] ) + t;
				LOG4CPP_TRACE( optLogger, "z-value of point " << i << ": " << fabs( p3dtrans( 2 ) ) );
				if ( fabs( p3dtrans( 2 ) ) > 1e-2 ) {
					LOG4CPP_TRACE( optLogger, "Points are NOT very coplanar" );
					//TODO ### UBITRACK_THROW( "Pose estimation requires four coplanar points" );
				}
				p3dAs2d.push_back( Math::Vector< 2 >( p3dtrans( 0 ), p3dtrans( 1 ) ) );
			}

			// compute homography
			Math::Matrix< 3, 3 > H;
			if ( p2d.size() > 4 )
				// copy first four elements to new vector
				H =  Calibration::homographyDLT( p3dAs2d, std::vector< Math::Vector< 2 > >( p2d.begin(), p2d.begin() + 4 ) );
			else
				H =  Calibration::homographyDLT( p3dAs2d, p2d );

			// compute initial pose from homography
			pose = Calibration::poseFromHomography( H, invK ) * Math::Pose( Math::Quaternion( P ), t );
			
			LOG4CPP_TRACE( optLogger, "Pose from homography (rotated): " << pose );
			Math::Matrix< 3, 3 > rotMat;
			pose.rotation().toMatrix( rotMat );
			LOG4CPP_TRACE( optLogger, "Rotation matrix (rotated): " << rotMat );
		}
	}
	
	// non-linear minimization
	Math::Matrix< 6, 6, double > covMatrix;
	if ( optimize )
	{
		residual = Calibration::optimizePose( pose, p2d, p3d, cam );
		LOG4CPP_DEBUG( optLogger, "Refined pose: " << pose << ", residual of 2D image measurements: " << residual);
	}
	else
	{
		residual = reprojectionError( p2d, p3d, pose, cam );
		LOG4CPP_DEBUG( optLogger, "NOT refined pose: " << pose << ", residual of 2D image measurements: " << residual);	
	}
	
	covMatrix = Calibration::singleCameraPoseError( pose, p3d, cam, residual );	
	residual = sqrt( residual / ( p2d.size() * 2) );
	
	return Math::ErrorPose( pose, covMatrix );
}



#endif // HAVE_LAPACK

} } // namespace Ubitrack::Calibration

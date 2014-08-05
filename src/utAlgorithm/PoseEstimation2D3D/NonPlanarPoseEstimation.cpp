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
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#include "NonPlanarPoseEstimation.h"
#include "../PoseEstimation3D3D/AbsoluteOrientation.h" // -> orientation estimation
#include <utMath/Blas1.h> // inner_product
#include <utMath/Blas2.h> // outer_product
#include <utMath/MatrixOperations.h> // matrix_inverse
#include <utMath/Geometry/PointProjection.h>
#include <utMath/Geometry/PointTransformation.h>

#include <vector>
#include <limits> // limit-max
#include <numeric> // std::accumulate
#include <algorithm> //std::transform

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <log4cpp/Category.hh>
static log4cpp::Category& optLogger( log4cpp::Category::getInstance( "Ubitrack.Calibration.2D3DPoseEstimation" ) );

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#endif
// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation2D3D {

#ifdef HAVE_LAPACK


namespace { // anonymous

/** @internal centroid calculation */
template< typename T >
Math::Vector< T, 3 > calculateCentroid( const std::vector< Math::Vector< T, 3 > > &points )
{
	Math::Vector< T, 3 > centroid = std::accumulate( points.begin(), points.end(), Math::Vector< T, 3 > ( 0, 0, 0 ) );
	return ( centroid / static_cast< T > ( points.size() ) );
}

/** @internal point shifting */
template< typename T >
Math::Vector< T, 3 > shiftToCenter( std::vector < Math::Vector< T, 3 > >& points )
{
	Math::Vector< T, 3 > centroid = calculateCentroid< T >( points );
	std::transform( points.begin(), points.end(), points.begin()
		, std::bind2nd( std::minus< Math::Vector< T, 3 > >( ), centroid ) );
	return Math::Vector< T, 3 >( centroid );
}

/** @internal calculates t-factor matrix */
template< typename T >
Math::Matrix< T, 3, 3 > calculateTFactorMatrix( const std::vector< Math::Matrix< T, 3, 3 > > &los )
{
	const std::size_t n = los.size();
	Math::Matrix< T, 3, 3 > tFactorMatrix = std::accumulate( los.begin() , los.end(), Math::Matrix< T, 3, 3 >::zeros() );
	tFactorMatrix /= n;
	tFactorMatrix = Math::Matrix< T, 3, 3 >::identity() - tFactorMatrix;
	tFactorMatrix = Math::invert_matrix( tFactorMatrix );
	return ( tFactorMatrix / n );
}

// this function got completely useless. the "old" absolute orientation is called instead
// /** @internal calculates the absolute orientation, old version was buggy @todo interface of absolute orientation should be changed.  */
// template< typename T >
// Math::Matrix< T, 3, 3 > absoluteOrientation( const std::vector< Math::Vector< T, 3 > > &pointsA, const std::vector< Math::Vector< T, 3 > >& pointsB )
// {
	// // Math::Pose pose = calculateAbsoluteOrientation( pointsA.begin(), pointsA.end(), pointsB.begin(), pointsB.end() );
	// Math::Pose pose = calculateAbsoluteOrientation( pointsA, pointsB );
	// return Math::Matrix< T, 3, 3 >( pose.rotation() );
	// // std::vector< Math::Matrix< T, 3, 3 > > matrices;
	// // matrices.reserve( pointsA.size() );
	// // std::transform( pointsA.begin(), pointsA.end(), pointsB.begin(), std::back_inserter( matrices )
		// // , Math::OuterProduct< 3, T >() );
	
	// // Math::Matrix< T, 3, 3 > B( Math::Matrix< T, 3, 3 >::zeros() );
	// // B = std::accumulate( matrices.begin(), matrices.end(), B );
	// // B /= pointsA.size();
	
	// // if ( Math::determinant( B ) < 0 )
		// // B *= -1;
	// // Math::Matrix< T, 3, 3 > U;
	// // Math::Vector< T, 3 > s;
	// // Math::Matrix< T, 3, 3 > Vt;
	
	// // boost::numeric::bindings::lapack::gesvd( 'A', 'A', B, s, U, Vt );
	// // return ublas::trans( ublas::prod(  U, Vt ) );
// }

/** @internal calculates the newest translation  */
template< typename T >
Math::Vector< T, 3 > estimateTranslation( const std::vector< Math::Matrix< T, 3, 3 > >& Vi, const Math::Matrix< T, 3, 3 > &Rot, const std::vector< Math::Vector< T, 3 > > &pointsObj, const Math::Matrix< T, 3, 3 > &TMatrix )
{
	std::vector< Math::Vector< T, 3 > > vec_tmp;
	vec_tmp.reserve( Vi.size() );
	Math::Geometry::transform_points( Rot, pointsObj.begin(), pointsObj.end(), std::back_inserter( vec_tmp ) );
	std::transform( Vi.begin(), Vi.end(), vec_tmp.begin(), vec_tmp.begin(), Math::Geometry::TransformPoint() );
	
	Math::Vector< T, 3 > translation = std::accumulate( vec_tmp.begin(), vec_tmp.end(), Math::Vector< T, 3 >( 0, 0, 0 ) );
	return ublas::prod( TMatrix, translation );
}

/** @internal calculation of object-space-error for std::transform */
template< typename T >
struct object_space_error
{
public:
	T operator() ( const Math::Matrix< T, 3, 3 >& matrix, const Math::Vector< T, 3 >& vec ) const
    {
		Math::Matrix< T, 3, 3 > tmp( Math::Matrix< T, 3, 3 >::identity() - matrix );
		Math::Vector< T, 3 > vec_tmp = ublas::prod( tmp, vec );
		return Math::inner_product( vec_tmp, vec_tmp );
	}
};

/** @internal calculation of summarized object space error */
template< typename T >
T calculateObjectSpaceError( const std::vector< Math::Matrix< T, 3, 3 > > &Vi, const std::vector< Math::Vector< T, 3 > > &points  )
{
	std::vector< T > ose;
	ose.reserve( Vi.size() );
	std::transform( Vi.begin(), Vi.end(), points.begin(), std::back_inserter( ose ), object_space_error< T >() );
	T sum( 0 );
	return std::accumulate( ose.begin(), ose.end(), sum );
}

/** @internal calculation of line-ofsight-projection-matrix, can be easily used with std::transform */
template< typename T >
struct lineOfSightProjectionMatrix
{
public:
    Math::Matrix< T, 3, 3 > operator() ( const Math::Vector< T, 3 > &vec1 ) const
    {
		const T d = Math::inner_product( vec1, vec1 );
		const Math::Vector< T, 3 > vec2 = vec1 * ( 1./ d );
		return Math::outer_product( vec1, vec2 );
    }
};

/** @internal */
template< typename T > 
bool estimatePose2D3D_impl( const std::vector< Math::Vector< T, 2 > >& p2D_in, Math::Pose& p, const std::vector< Math::Vector< T, 3 > >& p3D_in,  std::size_t &max_iter, T &max_error  )
{
	assert( max_iter > 0 );
	
	// algorithm expects homogeneous image coordinates ( x, y, 1 )^T
	// points are copied since they are altered during estimation
	std::vector< Math::Vector< T, 3 > > p2Dh;
	p2Dh.reserve( p2D_in.size() );
	// p2Dh.assign( p2D.begin(), p2D.end() ); // better: -> add the final coordinate
	typename std::vector< Math::Vector< T, 2 > >::const_iterator it = p2D_in.begin();
	typename std::vector< Math::Vector< T, 2 > >::const_iterator itEnd = p2D_in.end();
	for( ; it != itEnd; ++it )
		p2Dh.push_back( Math::Vector< T, 3 >( (*it)[ 0 ], (*it)[ 1 ], 1 ) );
	
	// translate the 3D object points to coordinate center
	// this could also be done outside. In repeated usage of this algorithm here is a minimal optimization possibility
	std::vector< Math::Vector< T, 3 > > p3D;
	p3D.reserve( p3D.size() );
	p3D.assign( p3D_in.begin(), p3D_in.end() );
	Math::Vector< T, 3 > center = shiftToCenter( p3D );
	

	// Line-of-sight projection matrices for each image point:
	// V_i = (p_i) * (p_i)^t / ((p_i)^t *(p_i))
	std::vector< Math::Matrix< T, 3, 3 > > Vi;
	Vi.reserve( p2Dh.size() );
	std::transform( p2Dh.begin(), p2Dh.end(), std::back_inserter( Vi ), lineOfSightProjectionMatrix< T >() );
	
	
	//Matrix as a factor for estimation of T
	const Math::Matrix< T, 3, 3 > TfactorMatrix ( calculateTFactorMatrix( Vi ) );

	T error_old = std::numeric_limits< T >::max();
	std::size_t iterations = 0;
	bool converged = false;
	for( ; iterations < max_iter; ++iterations )
	{
		// project image points into object space
		std::transform( Vi.begin(), Vi.end(), p2Dh.begin(),  p2Dh.begin(), Math::Geometry::TransformPoint() );

		// compute the optimal estimate of R
		shiftToCenter( p2Dh );
		
		// Object Points are already shifted at the beginning and never being changed
		Math::Matrix< T, 3, 3 > R;
		if( !PoseEstimation3D3D::estimateRotation_3D3D( p2Dh, R,  p3D ) )
			break;

		// compute new approximation of T
		Math::Vector< T, 3 > Tr = estimateTranslation( Vi, R, p3D, TfactorMatrix );
		
		//project objects points into camera coordinates
		Math::Matrix< T, 3, 4 > projection;
		ublas::subrange( projection, 0, 3, 0, 3 ) = R;
		ublas::column( projection, 3 ) = Tr;
		//apply the projection as transformation
		Math::Geometry::transform_points( projection, p3D.begin(), p3D.end(), p2Dh.begin() );
		
		T error_new = calculateObjectSpaceError( Vi, p2Dh );
		
		//check termination criteria 
		// converged = ( iterations >= max_iter ) || (error_new <  max_error ) ;
		converged = ( iterations >= max_iter ) || ( std::fabs( error_old - error_new ) <= max_error  );
		if( converged )
		{
			center = ublas::prod( R , center );
			Tr -= center;
			p = Math::Pose( Math::Quaternion( R ).normalize(), Tr  );
			error_old = error_new;
			break;
		}
		error_old = error_new;
		LOG4CPP_TRACE( optLogger, "object-space error " << error_new << " after " << iterations << " iterations." );
	}

	//check if algorithm terminated due to minimum error
	converged = ( iterations < max_iter ) && converged;
	
	// set return values
	max_error = error_old;
	max_iter = iterations;
	
	return converged;
}

} // anonymous-namespace

bool estimatePose6D_2D3D( const std::vector< Math::Vector2d >& p2D, Math::Pose& p, 
	const std::vector< Math::Vector3d >& p3D, std::size_t &max_iter, double &error )
{
	LOG4CPP_DEBUG( optLogger, "starting 2D-3D pose estimate with double values." );
	return estimatePose2D3D_impl( p2D, p, p3D, max_iter, error );
}

bool estimatePose6D_2D3D( const std::vector< Math::Vector2f >& p2D, Math::Pose& p, 
	const std::vector< Math::Vector3f >& p3D, std::size_t &max_iter, float &error )
{
	LOG4CPP_DEBUG( optLogger, "starting 2D-3D pose estimate with float values." );
	return estimatePose2D3D_impl( p2D, p, p3D, max_iter, error );
}

#endif // HAVE_LAPACK

} } } // namespace Ubitrack::Algorithm::PoseEstimation2D3D

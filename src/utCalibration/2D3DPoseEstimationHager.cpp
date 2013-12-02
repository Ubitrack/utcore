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

#include "2D3DPoseEstimationHager.h" 
 
#include <vector>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <functional>

#include <utMath/MatrixOperations.h>
#include <utMath/Functors/VectorNFunctors.h>
#include <utMath/Functors/Vector3Functors.h>

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <log4cpp/Category.hh>
static log4cpp::Category& optLogger( log4cpp::Category::getInstance( "Ubitrack.Calibration.2D3DPoseEstimation" ) );

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#endif
// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;

namespace Ubitrack { namespace Calibration {

#ifdef HAVE_LAPACK

template< typename T >
Math::Vector< 3, T > calculateCentroid( const std::vector< Math::Vector< 3, T > > &points )
{
	Math::Vector< 3, T > centroid = std::accumulate( points.begin(), points.end(), Math::Vector< 3, T > ( 0, 0, 0 ) );
	return ( centroid / static_cast< T > ( points.size() ) );
}

template< typename T >
Math::Vector< 3, T > shiftToCenter( std::vector < Math::Vector< 3, T > > &points )
{
	Math::Vector< 3, T > centroid = calculateCentroid< T >( points );
	std::transform( points.begin(), points.end(), points.begin()
		, std::bind2nd( Math::Functors::difference_vector< 3, T >(), centroid ) );
	return centroid;
}

template< typename T >
Math::Matrix< 3, 3, T > calculateTFactorMatrix( const std::vector< Math::Matrix< 3, 3, T > > &los )
{
	Math::Matrix< 3, 3, T > tFactorMatrix( Math::Matrix< 3, 3, T >::zeros() );
	tFactorMatrix = std::accumulate( los.begin() , los.end(), tFactorMatrix );
	tFactorMatrix /= los.size();
	tFactorMatrix = Math::Matrix< 3, 3, T >::identity() - tFactorMatrix;
	tFactorMatrix = Math::invert_matrix( tFactorMatrix );
	return ( tFactorMatrix / los.size() );
}

template< typename T >
Math::Matrix< 3, 3, T > absoluteOrientation( const std::vector< Math::Vector< 3, T > > &pointsA, const std::vector< Math::Vector< 3, T > > pointsB )
{

	std::vector< Math::Matrix< 3, 3, T > > matrices;
	matrices.reserve( pointsA.size() );
	std::transform( pointsA.begin(), pointsA.end(), pointsB.begin(), std::back_inserter( matrices )
		, Math::Functors::distinct_outer_product< 3, T >() );
	
	Math::Matrix< 3, 3, T > B( Math::Matrix< 3, 3, T >::zeros() );
	B = std::accumulate( matrices.begin(), matrices.end(), B );
	B /= pointsA.size();
	
	if ( Math::determinant( B ) < 0 )
		B *= -1;
	Math::Matrix< 3, 3, T > U;
	Math::Vector< 3, T > s;
	Math::Matrix< 3, 3, T > Vt;
	
	boost::numeric::bindings::lapack::gesvd( 'A', 'A', B, s, U, Vt );
	Math::Matrix< 3, 3, T > rotation ( ublas::trans( ublas::prod(  U, Vt ) ) );
	
	return rotation;
}

template< typename T >
void projectPoints( std::vector< Math::Vector< 3, T > > &pointsImg, Math::Matrix< 3, 3, T > &Rot, Math::Vector< 3, T > &Tr, const std::vector< Math::Vector< 3, T > > &pointsObj )
{
	//prepare Matrix
	Math::Matrix< 3, 4, T > projection;
	ublas::subrange( projection, 0, 3, 0, 3 ) = Rot;
	ublas::column( projection, 3 ) = Tr;
	
	//apply transformation
	std::transform( pointsObj.begin(), pointsObj.end(), pointsImg.begin()
		, std::bind1st( Math::Functors::transform3x4_vector3< T >(), projection ) );
}

template< typename T >
Math::Vector< 3, T > estimateTranslation( const std::vector< Math::Matrix< 3, 3, T > > &Vi, const Math::Matrix< 3, 3, T > &Rot, const std::vector< Math::Vector< 3, T > > &pointsObj, const Math::Matrix< 3, 3, T > &TMatrix )
{
	std::vector< Math::Vector< 3, T > > vec_tmp;
	vec_tmp.reserve( Vi.size() );
	std::transform( pointsObj.begin(), pointsObj.end(), std::back_inserter( vec_tmp )
		, std::bind1st( Math::Functors::transform3x3_vector3< T >(), Rot ) );
	
	std::transform( Vi.begin(), Vi.end(), vec_tmp.begin(), vec_tmp.begin()
		, Math::Functors::transform3x3_vector3< T >() );
	
	Math::Vector< 3, T > translation = std::accumulate( vec_tmp.begin(), vec_tmp.end(), Math::Vector< 3, T >( 0, 0, 0 ) );
	translation = ublas::prod( TMatrix, translation );
	return translation;
}

/** calculation of object-space-error for STL container*/
template< typename T >
struct object_space_error
	: public std::binary_function< Math::Matrix< 3, 3, T >, Math::Vector< 3, T >, T >
{
public:
	T operator() ( const Math::Matrix< 3, 3, T > &matrix, const Math::Vector< 3, T >& vec ) const
    {
		Math::Matrix< 3, 3, T > tmp( Math::Matrix< 3, 3, T >::identity() - matrix );
		Math::Vector< 3, T > vec_tmp = ublas::prod( tmp, vec );
		return ublas::inner_prod( vec_tmp, vec_tmp );
	}
};

template< typename T >
T calculateObjectSpaceError( std::vector< Math::Matrix< 3, 3, T > > &Vi, std::vector< Math::Vector< 3, T > > &points  )
{
	std::vector< T > ose;
	ose.reserve( Vi.size() );
	std::transform( Vi.begin(), Vi.end(), points.begin(), std::back_inserter( ose ), object_space_error< T >() );
	T sum( 0 );
	return std::accumulate( ose.begin(), ose.end(), sum );
}

template< typename T >
T abskernel( const std::vector< Math::Vector< 3, T > > &pointsObj, std::vector< Math::Vector< 3, T > > &pointsImg, std::vector< Math::Matrix< 3, 3, T > > &Vi, Math::Matrix< 3, 3, T > &Tmatrix, Math::Matrix< 3, 3, T > &Rot, Math::Vector< 3, T > &Tr )
{		

	std::transform( Vi.begin(), Vi.end(), pointsImg.begin(),  pointsImg.begin(), Math::Functors::transform3x4_vector3< T >() );

	// compute the optimal estimate of R
	shiftToCenter( pointsImg );
	// Object Points are already shifted at the beginning and never being changed
	Rot = absoluteOrientation( pointsObj, pointsImg );
	
	
	// compute new approximation of T
	Tr = estimateTranslation( Vi, Rot, pointsObj, Tmatrix );
	
	//project objects points into camera coordinates
	projectPoints( pointsImg, Rot, Tr,  pointsObj );
	

	T error = calculateObjectSpaceError( Vi, pointsImg );
	return error;
}


/** calculation of line-ofsight-projection-matrix for STL container */
template< typename T >
struct lineOfSightProjectionMatrix
	: public std::binary_function< Math::Vector< 3, T >, T, Math::Matrix< 3, 3, T > >
{
public:
    Math::Matrix< 3, 3, T > operator() ( const Math::Vector< 3, T > &v, const T& d ) const
    {
		Math::Matrix< 3, 3, T > matrix = ublas::outer_prod( v, v );
		matrix /= d;
		return matrix;
    }
};


/** \internal */
template< typename T > 
bool estimatePoseImpl( Math::Pose& p, const std::vector< Math::Vector< 3, T > >& p2D, std::vector< Math::Vector< 3, T > >& p3D, 
	const Math::Matrix< 3, 3, T >& cam, unsigned &nIterations, T &termination_error  )
{
	// image points should already be moved to the center and homogenized
	// now just copy the image points, since they will be changed during the algorithm
	std::vector< Math::Vector< 3, T > > p2Dh;
	p2Dh.reserve( p2D.size() );
	p2Dh.assign( p2D.begin(), p2D.end() );
	

	// move the 3D object points to coordinate center
	// this could also be done outside. In repreated usage of this algorithm here is a minimal optimization possibility
	Math::Vector< 3, T > center = shiftToCenter( p3D );
	

	//calculate the inner product for each vector
	std::vector< T > dotProduct;
	dotProduct.reserve( p2Dh.size() ); 
	std::transform( p2Dh.begin(), p2Dh.end(), std::back_inserter( dotProduct ), Math::Functors::inner_product< 3, T > () );
	
	
	//Line-of-sight projection matrices
	std::vector< Math::Matrix< 3, 3, T > > Vi;
	Vi.reserve( p2Dh.size() );
	std::transform( p2Dh.begin(), p2Dh.end(), dotProduct.begin(), std::back_inserter( Vi ), lineOfSightProjectionMatrix< T >() );
	
	
	//Matrix as a factor for estimation of T
	Math::Matrix< 3, 3, T > TfactorMatrix ( calculateTFactorMatrix( Vi ) );

		
	//starting the algorithm loop to estimate R and T
	Math::Matrix< 3, 3, T > R;
	Math::Vector< 3, T > Tr;
	T error_old, error_new;
	unsigned iterations = 1;
	bool converging = true;
	
	error_new = abskernel( p3D, p2Dh, Vi, TfactorMatrix, R, Tr );
	LOG4CPP_TRACE( optLogger, "\nError " << error_new << " after " << iterations << " iterations." );

	while( converging )
	{
		++iterations;
		error_old = error_new;
		
		error_new = abskernel( p3D, p2Dh, Vi, TfactorMatrix, R, Tr );
		
		//check termination criterias 
		converging = ( iterations <= nIterations ) && ( termination_error < error_new ) ;
		LOG4CPP_TRACE( optLogger, "\nError " << error_new << " after " << iterations << " iterations." );
	}

	center = ublas::prod( R , center );
	Tr -= center;
	p = Math::Pose( Math::Quaternion( R ).normalize(), Tr  );
	
	
	// if z-translation is negative invert pose
	if( Tr( 2 ) > 0 )
		p = ~p;

	//check if algorithm terminated due to minimum error
	bool converged = ( iterations < nIterations ) && !converging;
	
	// set return values
	termination_error = error_new;
	nIterations = iterations;
	
	return  converged;
}


bool estimatePose( Math::Pose& p, std::vector< Math::Vector< 3, float > > p2D,
	std::vector< Math::Vector< 3, float > > p3D, const Math::Matrix< 3, 3, float >& cam,
	unsigned &nIterations, float &error )
{
	LOG4CPP_DEBUG( optLogger, "starting Pose Estimate with float values." );
	return estimatePoseImpl( p, p2D, p3D, cam, nIterations, error );
}

bool estimatePose( Math::Pose& p, std::vector< Math::Vector< 3, double > > p2D,
	std::vector< Math::Vector< 3, double > > p3D, const Math::Matrix< 3, 3, double >& cam,
	unsigned &nIterations, double &error )
{
	LOG4CPP_DEBUG( optLogger, "starting Pose Estimate with double values." );
	return estimatePoseImpl( p, p2D, p3D, cam, nIterations, error );
}

#endif // HAVE_LAPACK

} } // namespace Ubitrack::Calibration

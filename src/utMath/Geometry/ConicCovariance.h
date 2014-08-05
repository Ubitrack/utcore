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
 * Covariance estimation for least-square conic parameter estimation.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_MATH_GEOMETRY_CONIC_COVARIANCE_ESTIMATION_H_INCLUDED__
#define __UBITRACK_MATH_GEOMETRY_CONIC_COVARIANCE_ESTIMATION_H_INCLUDED__

// Ubitrack
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Blas1.h> // norm_2
#include <utMath/Blas2.h> // outer_product
#include <utMath/MatrixOperations.h>


namespace Ubitrack { namespace Math { namespace Geometry {


/**
 * This function calculates the covariance of least-square estimated conic parameters.
 * 
 * The function uses an approach introduced from Kanatani 2008 in his article
 * "Statistical optimization for geometric fitting: Theoretical accuracy bound
 * and high order error analysis" ( @cite kanatani2008statistical ) in paragraph 3.4.
 *
 * @verbatim
@article{kanatani2008statistical,
  title={Statistical optimization for geometric fitting: Theoretical accuracy bound and high order error analysis},
  author={Kanatani, Kenichi},
  journal={International Journal of Computer Vision},
  volume={80},
  number={2},
  pages={167--188},
  year={2008},
  publisher={Springer}
} @endverbatim
 * 
 * @tparam InputIterator type of iterator to container of data points (2d point parameters)
 * @tparam T numeric type used within this function, depends on given conic vector value_type (e.g. \c double or \c float )
 * @param iBegin iterator to first element of 2d-points in a container describing the conic's outline.
 * @param iEnd iterator to final element of 2d-points in a container describing the conic's outline.
 * @param conic describes a conic as a 6-vector 
 * @param covariance this matrix will be set to the covariance of the estimated conic parameters
 */
template< typename InputIterator, typename T >
bool estimateCovariance( const InputIterator iBegin, const InputIterator iEnd
	, const Math::Vector< T, 6 > &conic
	// , Math::Matrix< conic::value_type, 6, 6 > &covariance )
	, Math::Matrix< typename std::iterator_traits< InputIterator >::value_type::value_type, 6, 6 > &covariance )
{
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type;
	typedef typename vector_type::value_type value_type;

	
	const std::size_t n = std::distance( iBegin, iEnd );
	const value_type n2 = Math::norm_2( conic );
	const value_type a = conic [ 0 ] / n2;
	const value_type b = conic [ 1 ] / n2;
	const value_type c = conic [ 2 ] / n2;
	const value_type d = conic [ 3 ] / n2;
	const value_type e = conic [ 4 ] / n2;
	const value_type f = conic [ 5 ] / n2;
	
	 /// @todo should be estimated and not set hard coded
	// const value_type stdDev = 1;//0.1;
	const value_type stdDev = 0.0001;
	
	Math::Matrix< value_type, 6, 6 > invCov = Math::Matrix< value_type, 6, 6 >::zeros();
	
	for( InputIterator it( iBegin ); it != iEnd; ++it )
	{
		const value_type x = (*it)[ 0 ];
		const value_type y = (*it)[ 1 ];
		 // some check that maybe included to avoid unusable input
		// if ( x != x || y != y )
			// continue;
			
		const value_type x2 = x*x;
		const value_type y2 = y*y;
		Math::Vector< T, 6 > t;
		t[ 0 ] = x2;
		t[ 1 ] = 2*x*y;
		t[ 2 ] = y2;
		t[ 3 ] = 2*x;
		t[ 4 ] = 2*y;
		t[ 5 ] = 1;
	
		Math::Matrix< value_type, 6, 6 > tmp = Math::outer_product( t, t );
		
		const value_type t4 = b*stdDev*x*y*4.0;
		const value_type t0 = b*(b*stdDev*(y2*4.0+x2*4.0)+d*stdDev*y*4.0+e*stdDev*x*4.0+a*stdDev*x*y*4.0+c*stdDev*x*y*4.0)+d*(d*stdDev*4.0+a*stdDev*x*4.0+b*stdDev*y*4.0)+e*(e*stdDev*4.0+b*stdDev*x*4.0+c*stdDev*y*4.0)+a*(t4+a*stdDev*x2*4.0+d*stdDev*x*4.0)+c*(t4+c*stdDev*y2*4.0+e*stdDev*y*4.0);

		tmp /= t0;
		invCov += tmp;
	}
	
	covariance = Math::pseudoInvert_matrix( invCov );
	
	return true;
}

	
}}} // namespace Ubitrack::Math::Geometry


#endif //__UBITRACK_MATH_GEOMETRY_CONIC_COVARIANCE_ESTIMATION_H_INCLUDED__

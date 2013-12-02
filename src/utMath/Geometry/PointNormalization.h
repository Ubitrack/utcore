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
 * @ingroup math
 * @file
 *
 * Several functions that can be applied to points or container of points.
 *
 * The header includes functions for common operations on points.
 * The functions can easily be combined to STL-containers like 
 * std::vectors, std::list, etc, that contain points.
 * representing conics or explicit conics.
 *
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __POINT_NORMALIZATION_H__
#define __POINT_NORMALIZATION_H__

// Ubitrack
#include <utMath/Vector.h>
#include <utMath/Matrix.h>

namespace Ubitrack { namespace Math { namespace Geometry {

/**
 * Computes the normalization parameters required for
 * numerical optimization (e.g. DLT) of a set of points.
 *
 * Mainly to be used internally by homography and projection matrix estimation.
 * The points have to be transformed according to p'=(p-shift)/scale
 *
 * @tparam N dimension of the point (usually 2 or 3 for 2D or 3D points)
 * @tparam T type of point ( e.g. \c double or \c float )
 * @tparam ForwardIterator type of iterator pointing to begin and end of container
 * @param iBegin iterator to the first element in a container/array
 * @param iEnd iterator to the final position of a container/array (usually \c .end() )
 * @param shift returns the mean value of the vectors.
 * @param scale returns the non-isotropic extension of the vectors around the mean
 */
template< std::size_t N, typename T, typename ForwardIterator >
void estimateNormalizationParameters( const ForwardIterator iBegin, const ForwardIterator iEnd
	, Math::Vector< N, T >& shift, Math::Vector< N, T >& scale )
{
	const std::size_t n_pts = std::distance( iBegin, iEnd );
	// compute mean and mean of square
	shift = Math::Vector< N, T >::zeros();
	scale = Math::Vector< N, T >::zeros();
	
	for ( ForwardIterator it( iBegin ); it < iEnd; ++it )
	{
		shift = shift + *it;
		scale = scale + boost::numeric::ublas::element_prod( *it, *it );
	}
	shift *= static_cast< T >( 1 ) / n_pts;
	scale *= static_cast< T >( 1 ) / n_pts;
	
	// compute standard deviation
	for ( std::size_t i = 0; i < N; i++ )
		scale( i ) = std::sqrt( scale( i ) - ( shift( i ) * shift( i ) ) );
}

/**
 * Generates the normalization matrix that corresponds to 
 * the given scale and shift parameters.
 *
 * This matrix can be used to easily generate a matrix
 * that can be applied to a container ( e.g. \c std::vector )
 * of Math::Vectors via a transformation Functor.
 *
 * @tparam N dimension of the point (usually 2 or 3 for 2D or 3D points )
 * @tparam T type of point ( e.g. \c double or \c float )
 * @param shift returns the mean value of the vectors.
 * @param scale returns the non-isotropic extension of the vectors around the mean
 * @param modInverse if true, the inverse matrix is returned
 * @return the (N+1)x(N+1) normalization matrix, corresponding to the given shift and scale parameters
 */

template< std::size_t N, typename T >
inline Math::Matrix< N+1, N+1, T > generateNormalizationMatrix( const Math::Vector< N, T >& shift, const Math::Vector< N, T >& scale, const bool modInverse )
{
	// compute correction matrix
	Math::Matrix< N+1, N+1, T > modMatrix = Math::Matrix< N+1, N+1, T >::zeros();
	modMatrix( N, N ) = static_cast< T >( 1 ); //homogeneous matrix :)
	
	if ( modInverse )
		for ( std::size_t i = 0; i < N; i++ )
		{
			modMatrix( i, i ) = scale( i );
			modMatrix( i, N ) = shift( i );
		}
	else
		for ( std::size_t i = 0; i < N; i++ )
		{
			modMatrix( i, i ) = static_cast< T >( 1 ) / scale( i );
			modMatrix( i, N ) = -modMatrix( i, i ) * shift( i );
		}
		
	return modMatrix;
};

}}} // namespace Ubitrack::Math::Geometry

#endif // __POINT_NORMALIZATION_H__
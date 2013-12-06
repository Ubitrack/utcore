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
 * Functors for common operations on 2-vectors
 *
 * The Functors can easily be applied to containers like 
 * std::vector of Math::Vectors< 2, T > using std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

 
 
#ifndef __H__VECTOR2_FUNCTORS__
#define __H__VECTOR2_FUNCTORS__

// std
#include <numeric>
#include <algorithm>
#include <functional>

// Ubitrack
#include <utMath/Vector.h>
#include <utMath/Matrix.h>

namespace Ubitrack { namespace Math { namespace Functors {

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * transforms a Math::Vector< T, 2 > by a given 3x3 (rotation-)matrix
 */
 
template< typename T >
struct transform3x3_vector2
	: public std::binary_function< Math::Matrix< 3, 3, T >, Math::Vector< T, 2 >, Math::Vector< T, 3 > >
{
public:
	/**
	 * @ingroup math
	 * transforms a Math::Vector< T, 2 > by a given 3x3 matrix
	 *
	 * @param mat the 3x3 rotation matrix
	 * @param vec the vector to be transformed
	 * @return transformed vector
	 */
	Math::Vector< T, 3 > operator() ( const Math::Matrix< 3, 3, T > &mat, const Math::Vector< T, 2 >& vec ) const
    {
		const T e1 = mat( 0, 0 ) * vec( 0 ) + mat( 0, 1 ) * vec( 1 ) + mat( 0, 2 );
		const T e2 = mat( 1, 0 ) * vec( 0 ) + mat( 1, 1 ) * vec( 1 ) + mat( 1, 2 );
		const T e3 = mat( 2, 0 ) * vec( 0 ) + mat( 2, 1 ) * vec( 1 ) + mat( 2, 2 );
		return Math::Vector< T, 3 > ( e1, e2, e3 );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Projects a Math::Vector< T, 3 > with a 3x4 projection matrix.
 */
 
template< typename T >
struct project3x3_vector2
	: public std::binary_function< Math::Matrix< 3, 3, T >, Math::Vector< T, 2 >, Math::Vector< T, 2 > >
{
public:
	/**
	 * @ingroup math
	 * Projects the Math::Vector< T, 3 > with a 3x4 projection matrix.
	 *
	 * @param projection the 3x4 projection matrix
	 * @param vec the 3-vector
	 * @return projected 3-vector
	 */
	Math::Vector< T, 2 > operator() ( const Math::Matrix< 3, 3, T > &projection, const Math::Vector< T, 2 > &vec ) const
    {
		const T e1 = projection( 0, 0 ) * vec( 0 ) + projection( 0, 1 ) * vec( 1 ) + projection( 0, 2 );
		const T e2 = projection( 1, 0 ) * vec( 0 ) + projection( 1, 1 ) * vec( 1 ) + projection( 1, 2 );
		const T e3 = projection( 2, 0 ) * vec( 0 ) + projection( 2, 1 ) * vec( 1 ) + projection( 2, 2 );
		return Math::Vector< T, 2 > ( e1/e3, e2/e3 );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Normalizes the Math::vector< 2, T >
 */
 
template< typename T >
struct normalize_vector2
	: public std::unary_function< Math::Vector< T, 2 >, Math::Vector< T, 2 > >
{
public:
	/**
	 * @ingroup math
	 * Normalizes the Math::vector< 2, T >
	 *
	 * @param v the 2-vector
	 * @return normalized vector
	 */
	Math::Vector< T, 2 > operator() ( const Math::Vector< T, 2 >& v ) const
    {
		const T norm = std::sqrt( v( 0 )*v( 0 )+v( 1 )*v( 1 ) );
		return ( v / norm );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Comparison function for image vectors (e.g. pixel coordinates).
 */
 
template< typename T >
struct is_less_vector2
	: public std::binary_function< Math::Vector< T, 2 >, Math::Vector< T, 2 >, bool >
{
public:
	/**
	 * @ingroup math
	 * Compares two distinct Math::Vector< T, 2 > for sorting algrotihms.
	 * Comaprison starts along the first dimension and, if equal,
	 * along the second dimension. Attention: This is not respecting 
	 * euclidean distance in any way.
	 *
	 * @param v1 the 1st 2-vector
	 * @param v2 the 2nd 2-vector
	 * @return true if first vector is_smaller than second
	 */
	bool operator() ( const Math::Vector< T, 2 >& v1, const Math::Vector< T, 2 >& v2 ) const
    {
		// compare the x-values
		if ( v1( 0 ) < v2( 0 ) )
			return true;
		//otherwise: comare y-values if x-values are equal
		return ( v1( 0 ) == v2( 0 ) )? v1( 1 ) < v2( 1 ) : false;
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the euclidean distance of two Math::Vector< T, 2 >.
 */
 
template< typename T >
struct euclidean_distance_vector2
	: public std::binary_function< Math::Vector< T, 2 >, Math::Vector< T, 2 >, T >
{
public:
 	/**
	 * @ingroup math
	 * Calculates the euclidean distance of two Math::Vector< T, 2 >.
	 *
	 * @param v1 the first 2-vector
	 * @param v2 the second 2-vector
	 * @return euclidean distance between v1 and v2
	 */
	T operator()( const Math::Vector< T, 2 > &v1, const Math::Vector< T, 2 > &v2 ) const
	{
		const T x = v1( 0 ) - v2( 0 );
		const T y = v1( 1 ) - v2( 1 );
		return std::sqrt( x*x + y*y );
	}
};

} } } // namespace Ubitrack::Math::Functors

#endif  // __H__VECTOR2_FUNCTORS__

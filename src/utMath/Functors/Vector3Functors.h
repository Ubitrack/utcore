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
 * Functors for common operations on 3-vectors
 *
 * The Functors can easily be applied to containers like 
 * std::vector of Math::Vectors< 3, T > using std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

 
 
#ifndef __H__VECTOR3_FUNCTORS__
#define __H__VECTOR3_FUNCTORS__

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
 * transforms a Math::Vector< 3, T > by a given 3x3 transformation-matrix.
 */
 
template< typename T >
struct transform3x3_vector3
	: public std::binary_function< Math::Matrix< 3, 3, T >, Math::Vector< 3, T >, Math::Vector< 3, T > >
{
public:
	/**
	 * @ingroup math
	 * transforms a Math::Vector< 3, T > by a given 3x3 transformation-matrix.
	 *
	 * @param mat the 3x3 transformation matrix
	 * @param vec the vector to be transformed
	 * @return transformed vector
	 */
	Math::Vector< 3, T > operator() ( const Math::Matrix< 3, 3, T > &mat, const Math::Vector< 3, T >& vec ) const
    {
		const T e1 = mat( 0, 0 ) * vec( 0 ) + mat( 0, 1 ) * vec( 1 ) + mat( 0, 2 ) * vec( 2 );
		const T e2 = mat( 1, 0 ) * vec( 0 ) + mat( 1, 1 ) * vec( 1 ) + mat( 1, 2 ) * vec( 2 );
		const T e3 = mat( 2, 0 ) * vec( 0 ) + mat( 2, 1 ) * vec( 1 ) + mat( 2, 2 ) * vec( 2 );
		return Math::Vector< 3, T > ( e1, e2, e3 );
	}
};


/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * transforms a Math::Vector< 3, T > by a given 3x4 transformation-matrix
 */
 
template< typename T >
struct transform3x4_vector3
	: public std::binary_function< Math::Matrix< 3, 4, T >, Math::Vector< 3, T >, Math::Vector< 3, T > >
{
public:
	/**
	 * @ingroup math
	 * transforms a Math::Vector< 3, T > by a given 3x4 transformation-matrix.
	 *
	 * @param mat the 3x4 transformation matrix
	 * @param vec the vector to be transformed
	 * @return transformed vector
	 */
	Math::Vector< 3, T > operator() ( const Math::Matrix< 3, 4, T > &mat, const Math::Vector< 3, T >& vec ) const
    {
		const T e1 = mat( 0, 0 ) * vec( 0 ) + mat( 0, 1 ) * vec( 1 ) + mat( 0, 2 ) * vec( 2 ) + mat( 0, 3 );
		const T e2 = mat( 1, 0 ) * vec( 0 ) + mat( 1, 1 ) * vec( 1 ) + mat( 1, 2 ) * vec( 2 ) + mat( 1, 3 );
		const T e3 = mat( 2, 0 ) * vec( 0 ) + mat( 2, 1 ) * vec( 1 ) + mat( 2, 2 ) * vec( 2 ) + mat( 2, 3 );
		return Math::Vector< 3, T > ( e1, e2, e3 );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Projects a Math::Vector< 3, T > with a 3x4 projection-matrix.
 */
 
template< typename T >
struct project3x4_vector3
	: public std::binary_function< Math::Matrix< 3, 4, T >, Math::Vector< 3, T >, Math::Vector< 2, T > >
{
public:
	/**
	 * @ingroup math
	 * Projects the Math::Vector< 3, T > with a 3x4 projection-matrix.
	 *
	 * @param projection the 3x4 projection matrix
	 * @param vec the 3-vector
	 * @return projected 3-vector
	 */
	Math::Vector< 2, T > operator() ( const Math::Matrix< 3, 4, T > &projection, const Math::Vector< 3, T > &vec ) const
    {
		const T e1 = projection( 0, 0 ) * vec( 0 ) + projection( 0, 1 ) * vec( 1 ) + projection( 0, 2 ) * vec( 2 ) + projection( 0, 3 );
		const T e2 = projection( 1, 0 ) * vec( 0 ) + projection( 1, 1 ) * vec( 1 ) + projection( 1, 2 ) * vec( 2 ) + projection( 1, 3 );
		const T e3 = 1 / ( projection( 2, 0 ) * vec( 0 ) + projection( 2, 1 ) * vec( 1 ) + projection( 2, 2 ) * vec( 2 ) + projection( 2, 3 ) );
		return Math::Vector< 2, T > ( e1*e3, e2*e3 );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Projects a Math::Vector< 3, T > with a 4x4 projection-matrix.
 */
 
template< typename T >
struct project4x4_vector3
	: public std::binary_function< Math::Matrix< 4, 4, T >, Math::Vector< 3, T >, Math::Vector< 3, T > >
{
public:
	/**
	 * @ingroup math
	 * Projects the Math::Vector< 3, T > with a 4x4 projection-matrix(e.g. a 4x4-homography).
	 *
	 * @param projection the 4x4 projection matrix
	 * @param vec the 3-vector
	 * @return projected 3-vector
	 */
	Math::Vector< 3, T > operator() ( const Math::Matrix< 4, 4, T > &projection, const Math::Vector< 3, T > &vec ) const
    {
		const T e1 = projection( 0, 0 ) * vec( 0 ) + projection( 0, 1 ) * vec( 1 ) + projection( 0, 2 ) * vec( 2 ) + projection( 0, 3 );
		const T e2 = projection( 1, 0 ) * vec( 0 ) + projection( 1, 1 ) * vec( 1 ) + projection( 1, 2 ) * vec( 2 ) + projection( 1, 3 );
		const T e3 = projection( 2, 0 ) * vec( 0 ) + projection( 2, 1 ) * vec( 1 ) + projection( 2, 2 ) * vec( 2 ) + projection( 2, 3 );
		const T e4 = 1 / ( projection( 3, 0 ) * vec( 0 ) + projection( 3, 1 ) * vec( 1 ) + projection( 3, 2 ) * vec( 2 ) + projection( 3, 3 ) );
		return Math::Vector< 3, T > ( e1*e4, e2*e4, e3*e4 );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the norm of a Math::Vector< 3, T >.
 */
 
template< typename T >
struct norm_vector3
	: public std::unary_function< Math::Vector< 3, T >, T >
{
public:
	/**
	 * @ingroup math
	 * Calculates the norm of a Math::Vector< 3, T >.
	 *
	 * @param v the 3-vector
	 * @return norm of the input 3-vector
	 */
	T operator() ( const Math::Vector< 3, T >& v ) const
    {
		return std::sqrt( v( 0 )*v( 0 )+v( 1 )*v( 1 )+v( 2 )*v( 2 ) );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Normalizes the Math::Vector< 3, T >.
 */
 
template< typename T >
struct normalize_vector3
	: public std::unary_function< Math::Vector< 3, T >, Math::Vector< 3, T > >
{
public:
	/**
	 * @ingroup math
	 * Normalizes the Math::Vector< 3, T >.
	 *
	 * @param v the 3-vector
	 * @return normalized vector
	 */
	Math::Vector< 3, T > operator() ( const Math::Vector< 3, T >& v ) const
    {
		const T norm = std::sqrt( v( 0 )*v( 0 )+v( 1 )*v( 1 )+v( 2 )*v( 2 ) );
		return ( v / norm );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Transforms a Math::Vector< 2 , T > into a homogenious 3-vector.
 */
 
template< typename T >
struct homogenize_vector3
	: public std::unary_function< Math::Vector< 2, T >, Math::Vector< 3, T > >
{
public:
	/**
	 * @ingroup math
	 * Transforms a Math::Vector< 2 , T > into a homogenious 3-vector by
	 * adding one dimension, which is set to one.
	 *
	 * @param vec the 2-vector
	 * @return homogenious 3-vector
	 */
	Math::Vector< 3, T > operator()( const Math::Vector< 2, T > &vec ) const
	{
		return Math::Vector< 3, T >( vec( 0 ), vec( 1 ), 1 );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Transforms a Math::Vector< 3 , T > into a 2-vector.
 */
 
template< typename T >
struct dehomogenize_vector3
	: public std::unary_function< Math::Vector< 3, T >, Math::Vector< 2, T > >
{
public:
	/**
	 * @ingroup math
	 * Transforms a Math::Vector< 3 , T > into a 2-vector by
	 * removing the last dimension.
	 *
	 * @param vec the 3-vector
	 * @return 2-vector
	 */
	Math::Vector< 2, T > operator()( const Math::Vector< 3, T > &vec ) const
	{
		return Math::Vector< 2, T >( vec( 0 ), vec( 1 ) );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the euclidean distance of two Math::Vector< 3, T >.
 */
 
template< typename T >
struct euclidean_distance_vector3
	: public std::binary_function< Math::Vector< 3, T >, Math::Vector< 3, T >, T >
{
public:
	/**
	 * @ingroup math
	 * Calculates the euclidean distance of two Math::Vector< 3, T >.
	 *
	 * @param v1 the 1st 3-vector
	 * @param v2 the 2nd 3-vector
	 * @return euclidean distance between v1 and v2
	 */
	T operator()( const Math::Vector< 3, T > &v1, const Math::Vector< 3, T > &v2 ) const
	{
		const T x = v1( 0 ) - v2( 0 );
		const T y = v1( 1 ) - v2( 1 );
		const T z = v1( 2 ) - v2( 2 );
		return std::sqrt( x*x + y*y + z*z );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Constructs the antisymmetric-skew matrix from a Math::Vector< 3, T >.
 */
 
template< typename T >
struct skew_matrix
	: public std::unary_function< Math::Vector< 3, T >, Math::Matrix< 3, 3, T > >
{
public:
	/**
	 * @ingroup math
	 * Constructs the antisymmetric-skew matrix from a Math::Vector< 3, T >.
	 *
	 * If the resulting skew matrix is applied to a 3-vector the resulting vector
	 * is the cross-product of the two vectors.
	 * Therefore the 3x3 skew-matrix s of vector v statisfies s * v = 0.
	 *
	 * @param v the 3-vector
	 * @return antisymmetric-skew matrix of v
	 */
	Math::Matrix< 3, 3, T > operator() ( const Math::Vector< 3, T >& v ) const
    {
		namespace ublas = boost::numeric::ublas;
		
		ublas::matrix< T, ublas::column_major > skew( 3, 3 );
		ublas::row( skew, 0 ) = Math::Vector< 3, T >( 0, -v( 2 ), v( 1 ) );
		ublas::row( skew, 1 ) = Math::Vector< 3, T >( v( 2 ), 0, -v( 0 ) );
		ublas::row( skew, 2 ) = Math::Vector< 3, T >( -v( 1 ), v( 0 ), 0 );
		return Math::Matrix< 3, 3, T > ( skew );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Constructs a new Math::vector< 3, T > as the cross-product of two Math::Vector< 3, T >.
 */
 
template< typename T >
struct cross_product
	: public std::binary_function< Math::Vector< 3, T >, Math::Vector< 3, T >, Math::Vector< 3, T > >
{
public:
	/**
	 * @ingroup math
	 * Constructs a new Math::Vector< 3, T > as the cross-product of two Math::Vector< 3, T >.
	 *
	 * @param a the 1st 3-vector
	 * @param b the 2nd 3-vector
	 * @return cross-product of the two input vectors
	 */
	Math::Vector< 3, T > operator() ( const Math::Vector< 3, T >& a, const Math::Vector< 3, T > &b ) const
    {
		const T e1 = a( 1 ) * b( 2 ) - a( 2 ) * b( 1 );
		const T e2 = a( 2 ) * b( 0 ) - a( 0 ) * b( 2 );
		const T e3 = a( 0 ) * b( 1 ) - a( 1 ) * b( 0 );
		return Math::Vector< 3 > ( e1, e2, e3 );
	}
};

} } } // namespace Ubitrack::Math::Functors

#endif  // __H__VECTOR_FUNCTORS__

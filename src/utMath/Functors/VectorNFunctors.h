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
 * Functors for common operations on vectors
 *
 * The Functors can easily be applied to containers like 
 * std::vector of Math::Vectors< N, T > using std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

 
 
#ifndef __H__VECTORN_FUNCTORS__
#define __H__VECTORN_FUNCTORS__

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
 * Calculates the inner product of a given Math::Vector< N, T >
 */
 
template< unsigned N, typename T >
struct inner_product
	: public std::unary_function< Math::Vector< N, T >, T >
{
public:
	/**
	 * @ingroup math
	 * Calculates the inner product of a given Math::Vector< N, T >
	 *
	 * @param v the input vector
	 * @return inner product of the vector v
	 */
    T operator() ( const Math::Vector< N, T > &vec ) const
    {
		return static_cast< T > ( boost::numeric::ublas::inner_prod( vec, vec ) );
    }
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the outer product of a Math::Vector< N, T >.
 */

template< unsigned N, typename T >
struct outer_product
	: public std::unary_function< Math::Vector< N, T >, Math::Matrix< N, N, T > >
{
public:
	/**
	 * @ingroup math
	 * Calculates the outer product of a Math::Vector< N, T >.
	 *
	 * The outer product is given as op = v * v^t.
	 *
	 * @param v the vector
	 * @return outer product of vector v
	 */
	Math::Matrix< N, N, T > operator() ( const Math::Vector< N, T >& vec ) const
    {
		return boost::numeric::ublas::outer_prod( vec, vec );
	}
};


/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the outer product of two Math::Vector< N, T >.
 */

template< unsigned N, typename T >
struct distinct_outer_product
	: public std::binary_function< Math::Vector< N, T >, Math::Vector< N, T >, Math::Matrix< N, N, T > >
{
public:
	/**
	 * @ingroup math
	 * Calculates the outer product of two Math::Vector< N, T > u and v.
	 *
	 * The outer product is given as op = u * v^t.
	 *
	 * @param u the first vector
	 * @param v the second vector
	 * @return outer product of v and u
	 */
	Math::Matrix< N, N, T > operator() ( const Math::Vector< N, T >& v, const Math::Vector< N, T >& u ) const
    {
		return boost::numeric::ublas::outer_prod( u, v );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the weighted outer product of a Math::Vector< N, T >.
 */
 
template< unsigned N, typename T >
struct outer_product_weighted
	: public std::binary_function< Math::Vector< N, T >, T, Math::Matrix< N, N, T > >
{
public:
	/**
	 * @ingroup math
	 * Calculates the weighted outer product of a Math::Vector< N, T >.
	 *
	 * The weighted outer product is given as op = ( v * v^t ) * w.
	 *
	 * @param v the vector
	 * @param w weight of the vector v.
	 * @return outer product of weighted vector v
	 */
    Math::Matrix< N, N, T > operator() ( const Math::Vector< N, T > &v, const T w ) const
    {
		Math::Vector< N, T > v_tmp = v * w;
		return boost::numeric::ublas::outer_prod( v_tmp, v );
    }
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Multiplies a scalar to each element of a given matrix.
 */

template< unsigned N, typename T >
struct multiply_matrix
	: public std::binary_function< T, Math::Matrix< N, N, T >, Math::Matrix< N, N, T > >
{
public:
	/**
	 * @ingroup math
	 * Multiplies a scalar to each element of a given matrix.
	 *
	 * @param scalar the scalar value
	 * @param matrix 
	 * @return scaled matrix
	 */
	Math::Matrix< N, N, T > operator() ( const T& scalar, const Math::Matrix< N, N, T >& matrix ) const
    {
		return ( scalar * matrix );
	}
};

/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the difference of two given vectors.
 */

template< unsigned N, typename T >
struct difference_vector
	: public std::binary_function< Math::Vector< N, T >, Math::Vector< N, T >, Math::Vector< N, T > >
{
public:
	/**
	 * @ingroup math
	 * Calculates the difference of two given vectors.
	 *
	 * @param vec1 the 1st n-vector
	 * @param vec2 the 2nd n-vector
	 * @return difference vector
	 */
	Math::Vector< N, T > operator() ( const Math::Vector< N, T >& vec1, const Math::Vector< N, T >& vec2 ) const
    {
		return vec1 - vec2;
	}
};


/**
 * @ingroup math
 * Functor Class for STL algorithms.
 * Calculates the norm of the Math::Vector< N, T >.
 */
 
template< unsigned N, typename T >
struct norm_vector
	: public std::unary_function< Math::Vector< N, T >, T >
{
public:
	/**
	 * @ingroup math
	 * Calculates the norm of the Math::Vector< N, T >.
	 *
	 * @param vec the n-vector
	 * @return norm of the n-vector
	 */
	T operator() ( const Math::Vector< N, T >& vec ) const
    {
		return boost::numeric::ublas::norm_2( vec );
	}
};

} } } // namespace Ubitrack::Math::Functors

#endif  // __H__VECTOR_FUNCTORS__

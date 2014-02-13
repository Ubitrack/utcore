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
 * Functors and functions for common matrix-vector operations. The
 * functions inside this file represent some common BLAS level 2
 * operations.
 *
 * The Functors can easily be applied to stl-containers like
 * \c std::vector, \c std::list, \c std::set, etc. containing
 * vectors and matrices data-structures.
 *
 * A common example would be \c std::vector< Math::Matrix3x3d >
 * and \c std::vector< Math::Vector3d > which are altered via
 * \c std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __H__BLAS_LEVEL_2__
#define __H__BLAS_LEVEL_2__

// #include "vector_traits.h"

#include "Vector.h"
#include "Matrix.h"

namespace Ubitrack { namespace Math {

// template< typename Vec1, typename Vec2 >
// struct vector_matrix_traits
// {};

// template< typename T1, std::size_t N1, typename T2, std::size_t N2 >
// struct vector_matrix_traits< Math::Vector< T1, N1 >, Math::Vector< T2, N2 > >
// {
	// typedef typename std::size_t size_type;
	// typedef typename T1 value_type;
	// typedef typename Math::Matrix< T1, N1, N2 > matrix_type;
	// static const size_type size1 = N1;
	// static const size_type size2 = N2;
// };

/**
 * @ingroup math
 * Functor class to calculate the outer product 
 * of two vectors using a recursive implementation.
 */
struct OuterProduct	
{
public:

	// template< typename VecType1, typename VecType2 >
	// typename Math::vector_matrix_traits< VecType1, VecType2 >::matrix_type operator() ( const VecType1& vec1, const VecType2& vec2 ) const
	// {
		// UBITRACK_STATIC_ASSERT( ( Math::has_fixed_storage< VecType1 >::value ), FIRST_VECTOR_EXPECTS_FIXED_STORAGE_TYPE );
		// UBITRACK_STATIC_ASSERT( ( Math::has_fixed_storage< VecType2 >::value ), SECOND_VECTOR_EXPECTS_FIXED_STORAGE_TYPE );
		
		// static const typename Math::vector_traits< VecType1 >::size_type size1 = Math::vector_traits< VecType1 >::size;
		// static const typename Math::vector_traits< VecType2 >::size_type size2 = Math::vector_traits< VecType2 >::size;
		// typedef typename Math::vector_matrix_traits< VecType1, VecType2 >::matrix_type MatType;
		// typedef typename Math::vector_matrix_traits< VecType1, VecType2 >::value_type value_type;
		
		// value_type matData [ size1*size2 ];
		// // attention: our matrix is column major but the constructor for
		// // accepting an array expects row-major representation -> :(
		// // therefore I chose row-major variant to fill in values
			
		// // call for column-major representation:
		// // outer_product_impl< size1, size1, size2 >() ( vec1, vec2, matData );
		
		// // call for row-major representation:
		// outer_product_impl< size2, size2, size1 >() ( vec2, vec1, matData );
		// return MatType( matData );
	// }
	
		/**
	 * @ingroup math functor
	 * Calculates the outer product of two vectors.
	 *
	 * @tparam T built-in type of input vectors ( e.g. \c double or \c float )
	 * @tparam N1 length of 1st input vector
	 * @tparam N2 length of 2nd input vector
	 * @param vec1 the 1st input vector 
	 * @param vec2 the 2nd input vector 
	 * @return outer product of the two vectors as matrix
	 */
	template< typename T, std::size_t N1, std::size_t N2 >
	typename Math::Matrix< T, N1, N2 > operator() ( const Math::Vector< T, N1 >& vec1, const Math::Vector< T, N2 >& vec2 ) const
	{
		// attention: our matrix is column major but the constructor for
		// accepting an array expects row-major representation -> :(
		// therefore I chose row-major variant to fill in values
		
		T matData [ N1*N2 ];		
		// call for row-major representation:
		outer_product_impl< N2, N2, N1 >() ( vec2, vec1, matData );
		return Math::Matrix< T, N1, N2 >( matData );
	}
	
protected:
	/**
	 * @ingroup math functor
	 * Internal functor that implements the calculation of the
	 * outer product of two vectors.
	 *
	 * The implementation uses recursive functional programming,
	 * that optimizes calculations at compile time.
	 *
	 * @tparam M_NOW value that signs the actual element of the first input vector
	 * @tparam M value that defines the matrix rows
	 * @tparam N value that defines the matrix columns
	 */
	template< std::size_t M_NOW, std::size_t M, std::size_t N >
	struct outer_product_impl
	{
		template< typename VecType1, typename VecType2, typename MatType  >
		void operator() ( const VecType1& vec1, const VecType2& vec2, MatType mat[] ) const
		{	
			mat[ (M*(N-1)) + (M_NOW-1) ] = vec1[ M_NOW-1 ] * vec2[ N-1 ];
			outer_product_impl< M_NOW-1, M, N  >()( vec1, vec2, mat );
		}
	};
	
	/**
	 * Partial specialization of functor for setting the next column to compute.
	 *
	 * @tparam M value that defines the matrix rows
	 * @tparam N value that defines the matrix columns
	 */
	template< std::size_t M, std::size_t N  >
	struct outer_product_impl< 0, M, N >
	{
		template< typename VecType1, typename VecType2, typename MatType  >
		void operator() ( const VecType1& vec1, const VecType2& vec2, MatType mat[] ) const
		{	
			outer_product_impl< M, M, N-1 >()( vec1, vec2, mat );
		}
	};
	
	/**
	 * Partial specialization of functor that signs the final element
	 *
	 * @tparam M value that defines the matrix rows
	 * @tparam N value that defines the matrix columns
	 */
	template< std::size_t V1, std::size_t V2 >
	struct outer_product_impl< V1, V2, 0 >
	{
		template< typename VecType1, typename VecType2, typename MatType  >
		void operator() ( const VecType1&, const VecType2&, MatType mat[] ) const
		{	
			return;
		}
	};
};

} } // namespace Ubitrack::Math

#endif //__H__BLAS_LEVEL_2__

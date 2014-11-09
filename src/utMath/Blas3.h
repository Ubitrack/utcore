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
 * Functors and functions for common matrix-matrix operations. The
 * functions inside this file represent some common BLAS level 3
 * operations.
 *
 * The Functors can easily be applied to stl-containers like
 * \c std::vector, \c std::list, \c std::set, etc. containing
 * matrices data-structures.
 *
 * A common example would be \c std::vector< Math::Matrix2x3d >
 * and \c std::vector< Math::Matrix3x4d > which are altered via
 * \c std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __UBITRACK_MATH_BLAS_LEVEL_3_H__
#define __UBITRACK_MATH_BLAS_LEVEL_3_H__

#include "Util/matrix_traits.h"
#include "Geometry/container_traits.h"

#include "Matrix.h"

#include <iterator> // std::iterator_traits
#include <algorithm> // std::transform

namespace Ubitrack { namespace Math {


/**
 * @internal
 * @ingroup math
 * Functor template class to calculate a matrix-matrix product. Needs
 * partial specialization to define some action for certain types.
 */
#ifndef __UBITRACK_MATH_BLAS_LEVEL_2_H__
template< typename LeftType, typename RightType, typename ReturnType > struct Product;
#endif

/// @internal Partial specialization for a matrix-matrix product
template< typename M1T, std::size_t M1R, std::size_t M1C, typename M2T, std::size_t M2C, std::size_t M2R >
struct Product< Math::Matrix< M1T, M1R, M1C >, Math::Matrix< M2T, M2R, M2C >, Math::Matrix< M1T, M1R, M2C > >
{
public:

	template< typename MatrixTypeLeft, typename MatrixTypeRight, typename RetMatrixType >
	void operator() ( const MatrixTypeLeft& lhs, const MatrixTypeRight& rhs, RetMatrixType & result ) const
	{
		// UBITRACK_STATIC_ASSERT( ( Math::Util::matrix_traits< MatrixType >::storage_type == Util::fixed_storage_tag ), EXPECTS_MATRIX_OF_FIXED_STORAGE_TYPE );
		// UBITRACK_STATIC_ASSERT( ( Math::Util::vector_traits< VectorType >::storage_type == Util::fixed_storage_tag ), EXPECTS_VECTOR_OF_FIXED_STORAGE_TYPE );
		typedef typename Math::Util::matrix_traits< MatrixTypeLeft >::value_type value_type;
		
		//determine rows (size1) and columns (size2) of input/result matrix
		static const typename Math::Util::matrix_traits< MatrixTypeLeft >::size_type size1 = Math::Util::matrix_traits< MatrixTypeLeft >::size1;
		static const typename Math::Util::matrix_traits< MatrixTypeLeft >::size_type size2 = Math::Util::matrix_traits< MatrixTypeLeft >::size2;
		static const typename Math::Util::matrix_traits< MatrixTypeRight >::size_type size3 = Math::Util::matrix_traits< MatrixTypeRight >::size2;

		// call for column-major representation (there is no row_major call yet)
		mat_mat_product_impl_forward< size1, size2, size3, 0, 0, 0 >() ( lhs, rhs, result );
		//mat_mat_product_impl_backward< size1, size2, size1, size2 >() ( lhs, rhs, result );
	}
	
	template< typename MatrixTypeLeft, typename MatrixTypeRight >
	Math::Matrix< M1T, M1R, M2C > operator() ( const MatrixTypeLeft& lhs, const MatrixTypeRight& rhs ) const
	{
		Math::Matrix< M1T, M1R, M2C > result;
		this->operator()( lhs, rhs, result );
		return result;
	}
	
protected:
	/**
	 * @internal
	 * Internal functor that implements a matrix-matrix product.
	 *
	 * The implementation uses recursive functional programming,
	 * that optimizes calculations at compile time.
	 *
	 * @note it sums up elements starting from the first element.
	 *
	 * @tparam M1_MAX value that signs the first dimension of the 1st matrix (rows)
	 * @tparam N1_MAX value that signs the second dimension of the 1st matrix (columns)
	 * @tparam N2_MAX value that signs the second dimension of the 2nd matrix (columns, rows should be equal to the columns of the 1st matrix )
	 * @tparam M an integer that defines the row within the the resulting matrix
	 * @tparam N an integer that defines the column within the the resulting matrix
	 * @tparam STEP integer that defines the index to the input matrix elements which are processed for the m-n element of the resulting matrix.
	 */
	template< std::size_t M1_MAX, std::size_t N1_MAX, std::size_t N2_MAX, std::size_t M, std::size_t N, std::size_t STEP >
	struct mat_mat_product_impl_forward
	{
		template< typename MatTypeLeft, typename MatTypeRight, typename ReturnMatType  >
		void operator() ( const MatTypeLeft& mat1, const MatTypeRight& mat2, ReturnMatType &resultMat ) const
		{
			// std::cout << " Matrix " << M << " x " << N << " step " << STEP << " :";
			
			const typename Math::Util::matrix_traits< MatTypeLeft >::value_type value1 = Math::Util::matrix_traits< MatTypeLeft >::ptr( mat1 )[ (STEP*M1_MAX) + M ];
			const typename Math::Util::matrix_traits< MatTypeRight >::value_type value2 = Math::Util::matrix_traits< MatTypeRight >::ptr( mat2 )[ (N*N1_MAX) + STEP ];
			Math::Util::matrix_traits< ReturnMatType >::ptr( resultMat )[ (N*M1_MAX) + M ] += value1 * value2;
			// std::cout << value1 << " x  " << value2 << "\n";
			mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, M, N, STEP+1 >()( mat1, mat2, resultMat );
		}
	};
	
	/// @internal Partial specialization of functor that signs the first element (first row)
	template< std::size_t M1_MAX, std::size_t N1_MAX, std::size_t N2_MAX, std::size_t M, std::size_t N >//, std::size_t STEP >
	struct mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, M, N, 0 >
	{
		template< typename MatTypeLeft, typename MatTypeRight, typename ReturnMatType  >
		void operator() ( const MatTypeLeft& mat1, const MatTypeRight& mat2, ReturnMatType &resultMat ) const
		{
			
			// std::cout << " Matrix new init " << M << " x " << N << " : ";
			const typename Math::Util::matrix_traits< MatTypeLeft >::value_type value1 = Math::Util::matrix_traits< MatTypeLeft >::ptr( mat1 )[ M ];
			const typename Math::Util::matrix_traits< MatTypeRight >::value_type value2 = Math::Util::matrix_traits< MatTypeRight >::ptr( mat2 )[ (N*N1_MAX) ];
			Math::Util::matrix_traits< ReturnMatType >::ptr( resultMat )[ (N*M1_MAX) + M ] = value1 * value2;
			// std::cout << value1 << " x  " << value2 << "\n";
			mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, M, N, 1 >()( mat1, mat2, resultMat );
		}
	};
	
	
	/// @internal Partial specialization of functor that signs the next column, resets the row as well
	template< std::size_t M1_MAX, std::size_t N1_MAX, std::size_t N2_MAX, std::size_t N >
	struct mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, M1_MAX, N, 0 >
	{
		template< typename MatTypeLeft, typename MatTypeRight, typename ReturnMatType >
		void operator() ( const MatTypeLeft& mat1, const MatTypeRight& mat2, ReturnMatType &resultMat ) const
		{
			
			mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, 0, N+1 , 0 >()( mat1, mat2, resultMat );
		}
	};
	
	/// @internal Partial specialization of functor that signs the next row
	template< std::size_t M1_MAX, std::size_t N1_MAX, std::size_t N2_MAX, std::size_t M, std::size_t N >
	struct mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, M, N, N1_MAX >
	{
		template< typename MatTypeLeft, typename MatTypeRight, typename ReturnMatType  >
		void operator() ( const MatTypeLeft& mat1, const MatTypeRight& mat2, ReturnMatType &resultMat ) const
		{
			mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, M+1, N , 0 >()( mat1, mat2, resultMat );
		}
	};
	
	/// @internal Partial specialization of functor that signs the final element
	template< std::size_t M1_MAX, std::size_t N1_MAX, std::size_t N2_MAX >
	struct mat_mat_product_impl_forward< M1_MAX, N1_MAX, N2_MAX, 0, N2_MAX, 0 >
	{
		template< typename MatTypeLeft, typename MatTypeRight, typename ReturnMatType  >
		void operator() ( const MatTypeLeft& mat1, const MatTypeRight& mat2, ReturnMatType &resultMat ) const
		{
			return;
		}
	};
};

// include guard needed since same function is within BLAS2 header
#ifndef __UBITRACK_MATH_BLAS_LEVEL_2_H__
template< typename Type1, typename Type2, typename ReturnType >
inline void product( const Type1& lhs, const Type2& rhs, ReturnType &result )
{
	Product< Type1, Type2, ReturnType >()( lhs, rhs, result );
}

template< typename InputIterator1, typename InputIterator2, typename OutputIterator >
inline void product( const InputIterator1 iBegin1, const InputIterator1 iEnd1, const InputIterator2 iBegin2, OutputIterator iOut )
{
	typedef typename std::iterator_traits< InputIterator1 >::value_type input_type1;
	typedef typename std::iterator_traits< InputIterator2 >::value_type input_type2;
	typedef typename Ubitrack::Util::container_traits< OutputIterator >::value_type output_type;

	std::transform( iBegin1, iEnd1, iBegin2, iOut, Product< input_type1, input_type2, output_type > () );
}
#endif

} } // namespace Ubitrack::Math

#endif //__UBITRACK_MATH_BLAS_LEVEL_3_H__

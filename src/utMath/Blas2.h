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


#ifndef __UBITRACK_MATH_BLAS_LEVEL_2_H__
#define __UBITRACK_MATH_BLAS_LEVEL_2_H__

#include "Util/vector_traits.h"
#include "Util/matrix_traits.h"
#include "Geometry/container_traits.h"
#include "Stochastic/identity_iterator.h"

#include "Vector.h"
#include "Matrix.h"

#include <iterator> // std::iterator_traits
#include <algorithm> // std::transform

namespace Ubitrack { namespace Math {


/**
 * @internal
 * @ingroup math
 * Functor class template to calculate a outer product. Needs partial spevialization
 * to define action on certain input types.
 */
template< typename LeftType, typename RightType, typename ReturnType > struct OuterProduct;

 
/// @internal Partial sepcialization for a matrix-vector product
template< typename VT, typename std::size_t N1, std::size_t N2 >
struct OuterProduct< Math::Vector< VT, N1 >, Math::Vector< VT, N2 >, Math::Matrix< VT, N1, N2 > >
{
public:
	
	template< typename VectorType1, typename VectorType2, typename RetType >
	void operator() ( const VectorType1& vec1, const VectorType2& vec2, RetType & result ) const
	{
	
		//determine rows (size1) and columns (size2)
		static const typename Math::Util::vector_traits< VectorType1 >::size_type size1 = Math::Util::vector_traits< VectorType1 >::size;
		static const typename Math::Util::vector_traits< VectorType2 >::size_type size2 = Math::Util::vector_traits< VectorType2 >::size;

		// attention: our matrix is column major but the constructor for
		// accepting an array expects row-major representation -> :(
		// therefore I chose row-major variant to fill in values
		// changed again since new version uses raw array to fill in values...
		// outer_product_impl< size2, size2, size1 >() ( vec2, vec1, ptr );
		
		
		VT *ptr = Math::Util::matrix_traits< RetType >::ptr( result );
		
		// call for column-major representation:
		outer_product_impl< size1, size1, size2 >() ( vec1, vec2, ptr );
	}
	
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
	
	template< typename VectorType1, typename VectorType2 >
	typename Math::Matrix< VT, N1, N2 > operator() ( const VectorType1& vec1, const VectorType2& vec2 ) const
	{
	
		//determine rows (size1) and columns (size2)
		static const typename Math::Util::vector_traits< VectorType1 >::size_type size1 = Math::Util::vector_traits< VectorType1 >::size;
		static const typename Math::Util::vector_traits< VectorType2 >::size_type size2 = Math::Util::vector_traits< VectorType2 >::size;

		
		Math::Matrix< VT, size1, size2 > result;
		this->operator()( vec1, vec2, result );
		return result;
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
		template< typename VectorType1, typename VectorType2, typename MatType  >
		void operator() ( const VectorType1& vec1, const VectorType2& vec2, MatType mat[] ) const
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
		template< typename VectorType1, typename VectorType2, typename MatType  >
		void operator() ( const VectorType1& vec1, const VectorType2& vec2, MatType mat[] ) const
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
		template< typename VectorType1, typename VectorType2, typename MatType  >
		void operator() ( const VectorType1&, const VectorType2&, MatType mat[] ) const
		{	
			return;
		}
	};
};



/**
 * @internal
 * @ingroup math
 * Functor template class to calculate a product of matrix of vector. Needs
 * partial specialization to define some action for certain types.
 */
template< typename LeftType, typename RightType, typename ReturnType > struct Product;


/// @internal Partial sepcialization for a matrix-vector product
template< typename MT, typename std::size_t MC, std::size_t MR, typename VT, std::size_t VR >
struct Product< Math::Matrix< MT, MR, MC >, Math::Vector< VT, VR >, Math::Vector< VT, MR > >
{
public:

	template< typename MatrixType, typename VectorType, typename RetType >
	void operator() ( const MatrixType& lhs, const VectorType& rhs, RetType & result ) const
	{
		// UBITRACK_STATIC_ASSERT( ( Math::Util::matrix_traits< MatrixType >::storage_type == Util::fixed_storage_tag ), EXPECTS_MATRIX_OF_FIXED_STORAGE_TYPE );
		// UBITRACK_STATIC_ASSERT( ( Math::Util::vector_traits< VectorType >::storage_type == Util::fixed_storage_tag ), EXPECTS_VECTOR_OF_FIXED_STORAGE_TYPE );
		typedef typename Math::Util::matrix_traits< MatrixType >::value_type value_type;
		
		//determine rows (size1) and columns (size2)
		static const typename Math::Util::matrix_traits< MatrixType >::size_type size1 = Math::Util::matrix_traits< MatrixType >::size1;
		static const typename Math::Util::matrix_traits< MatrixType >::size_type size2 = Math::Util::matrix_traits< MatrixType >::size2;

		// call for column-major representation (there is no row_major call yet)
		// mat_vec_product_impl_forward< size1, size2, 0, 0 >() ( lhs, rhs, result );
		mat_vec_product_impl_backward< size1, size2, size1, size2 >() ( lhs, rhs, result );
	}
	
	template< typename MatrixType, typename VectorType >
	Math::Vector< VT, MR > operator() ( const MatrixType& lhs, const VectorType& rhs ) const
	{
		Math::Vector< VT, MR > result;
		this->operator()( lhs, rhs, result );
		return result;
	}
	
protected:
	/**
	 * @internal
	 * Internal functor that implements a matrix-vector product.
	 *
	 * The implementation uses recursive functional programming,
	 * that optimizes calculations at compile time.
	 *
	 * @note it sums up elements starting from the first element.
	 *
	 * @tparam M_MAX value that signs the first dimension of the matrix (rows)
	 * @tparam N_MAX value that signs the second dimension of the matrix (columns)
	 * @tparam M value that defines the current matrix/result vector row
	 * @tparam N value that defines the current matrix/vector column 
	 */
	template< std::size_t M_MAX, std::size_t N_MAX, std::size_t M, std::size_t N >
	struct mat_vec_product_impl_forward
	{
		template< typename MatType, typename VecType, typename ReturnType  >
		void operator() ( const MatType& mat, const VecType& vec1, ReturnType &vec2 ) const
		{
			vec2[ M ] += Math::Util::matrix_traits< MatType >::ptr( mat )[ (N*(M_MAX) ) + M ] * vec1[ N ];
			mat_vec_product_impl_forward< M_MAX, N_MAX, M, N+1 >()( mat, vec1, vec2 );
		}
	};
	
	/// @internal Partial specialization of functor that signs the first element (first row)
	template< std::size_t M_MAX, std::size_t N_MAX, std::size_t M >
	struct mat_vec_product_impl_forward< M_MAX, N_MAX, M, 0 >
	{
		template< typename MatType, typename VecType, typename ReturnType  >
		void operator() ( const MatType& mat, const VecType& vec1, ReturnType &vec2 ) const
		{
			vec2[ M ] = Math::Util::matrix_traits< MatType >::ptr( mat )[ M ] * vec1[ 0 ];
			mat_vec_product_impl_forward< M_MAX, N_MAX, M, 1 >()( mat, vec1, vec2 );
		}
	};
	
	
	/// @internal Partial specialization of functor that signs the next row
	template< std::size_t M_MAX, std::size_t N_MAX, std::size_t M >
	struct mat_vec_product_impl_forward< M_MAX, N_MAX, M, N_MAX >
	{
		template< typename MatType, typename VecType, typename ReturnType >
		void operator() ( const MatType& mat, const VecType& vec1, ReturnType &vec2 ) const
		{
			mat_vec_product_impl_forward< M_MAX, N_MAX, M+1, 0 >()( mat, vec1, vec2 );
		}
	};
	
	/// @internal Partial specialization of functor that signs the final element
	template< std::size_t M_MAX, std::size_t N_MAX >
	struct mat_vec_product_impl_forward< M_MAX, N_MAX, M_MAX, 0 >
	{
		template< typename MatType, typename VecType, typename ReturnType >
		void operator() ( const MatType& mat, const VecType&, const ReturnType& ) const
		{
			return;
		}
	};

	/**
	 * @internal
	 * Internal functor that implements a matrix-vector product.
	 *
	 * The implementation uses recursive functional programming,
	 * that optimizes calculations at compile time.
	 *
	 * @note it sums up elements starting from the final element.
	 *
	 * @tparam M_MAX value that signs the first dimension of the matrix (rows)
	 * @tparam N_MAX value that signs the second dimension of the matrix (columns)
	 * @tparam M value that defines the current matrix/result vector row
	 * @tparam N value that defines the current matrix/vector column 
	 */
	template< std::size_t M_MAX, std::size_t N_MAX, std::size_t M, std::size_t N >
	struct mat_vec_product_impl_backward
	{
		template< typename MatType, typename VecType, typename ReturnType  >
		void operator() ( const MatType& mat, const VecType& vec1, ReturnType &vec2 ) const
		{
			vec2[ M-1 ] += Math::Util::matrix_traits< MatType >::ptr( mat )[ ((N-1)*(M_MAX) ) + (M-1) ] * vec1[ N-1 ];
			mat_vec_product_impl_backward< M_MAX, N_MAX, M, N-1 >()( mat, vec1, vec2 );
		}
	};
	
	/// @internal Partial specialization of functor that signs the first element (final row)
	template< std::size_t M_MAX, std::size_t N_MAX, std::size_t M >
	struct mat_vec_product_impl_backward< M_MAX, N_MAX, M, N_MAX >
	{
		template< typename MatType, typename VecType, typename ReturnType  >
		void operator() ( const MatType& mat, const VecType& vec1, ReturnType &vec2 ) const
		{
			vec2[ M-1 ] = Math::Util::matrix_traits< MatType >::ptr( mat )[ ((N_MAX-1)*(M_MAX) ) + (M-1) ] * vec1[ N_MAX-1 ];
			mat_vec_product_impl_backward< M_MAX, N_MAX, M, N_MAX-1 >()( mat, vec1, vec2 );
		}
	};
	
	/// @internal Partial specialization of functor that signs the next row
	template< std::size_t M_MAX, std::size_t N_MAX, std::size_t M >
	struct mat_vec_product_impl_backward< M_MAX, N_MAX, M, 0 >
	{
		template< typename MatType, typename VecType, typename ReturnType >
		void operator() ( const MatType& mat, const VecType& vec1, ReturnType &vec2 ) const
		{
			mat_vec_product_impl_backward< M_MAX, N_MAX, M-1, N_MAX >()( mat, vec1, vec2 );
		}
	};
	
	/// @internal Partial specialization of functor that signs the final element (end of first row reached)
	template< std::size_t M_MAX, std::size_t N_MAX >
	struct mat_vec_product_impl_backward< M_MAX, N_MAX, 0, N_MAX >
	{
		template< typename MatType, typename VecType, typename ReturnType >
		void operator() ( const MatType& mat, const VecType&, const ReturnType& ) const
		{
			return;
		}
	};
};

/**
 * @ingroup math
 * @brief A function that calculates the outer product of two vectors.
 * 
 * This function calculates the outer product of two vectors \f$ u \f$
 * and \f$ v \f$ , each one consisting of n elements, as
 * \f$ u \cdot v^T  \f$ .
 *
 * This function template wraps a call to the \b OuterProduct
 * \b functor. The input vectors for calculating the outer product can be
 * of any dimension. 
 *
 * @tparam VecType1 the type of the 1st input vectors.
 * @tparam VecType1 the type of the 2nd input vectors.
 * @param lhs the \b 1st input vector
 * @param rhs the \b 2nd input vector
 * @return the outer product of the two vectors as a Matrix
 */
template< typename VectorType1, typename VectorType2 >
inline Math::Matrix< typename Math::Util::vector_traits< VectorType1 >::value_type, Math::Util::vector_traits< VectorType1 >::size, Math::Util::vector_traits< VectorType2 >::size >
outer_product( const VectorType1& lhs, const VectorType2& rhs )
{
	typedef typename Math::Util::vector_traits< VectorType1 >::value_type value_type;
	static const typename Math::Util::vector_traits< VectorType1 >::size_type size1 = Math::Util::vector_traits< VectorType1 >::size;
	static const typename Math::Util::vector_traits< VectorType2 >::size_type size2 = Math::Util::vector_traits< VectorType2 >::size;
		
	typedef typename Math::Matrix< value_type, size1, size2 > return_type ;
	
	return OuterProduct< VectorType1, VectorType2, return_type >()( lhs, rhs );
}

template< typename VectorType1, typename VectorType2, typename ReturnType >
inline void outer_product( const VectorType1& lhs, const VectorType2& rhs, ReturnType &result )
{
	OuterProduct< VectorType1, VectorType2, ReturnType >()( lhs, rhs, result );
}


template< typename InputIterator1, typename InputIterator2, typename OutputIterator >
inline void outer_product( const InputIterator1 iBegin1, const InputIterator1 iEnd1, const InputIterator2 iBegin2, OutputIterator iOut )
{
	typedef typename std::iterator_traits< InputIterator1 >::value_type input_type1;
	typedef typename std::iterator_traits< InputIterator2 >::value_type input_type2;
	typedef typename Ubitrack::Util::container_traits< OutputIterator >::value_type output_type;

	std::transform( iBegin1, iEnd1, iBegin2, iOut, OuterProduct< input_type1, input_type2, output_type > () );
}

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

} } // namespace Ubitrack::Math

#endif //__UBITRACK_MATH_BLAS_LEVEL_2_H__

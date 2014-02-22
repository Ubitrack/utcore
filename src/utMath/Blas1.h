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
 * Functors and functions for common operations on vectors. The
 * functions inside this file represent common BLAS level 1
 * operations.
 *
 * The Functors can easily be applied to stl-containers like
 * \c std::vector, \c std::list, \c std::set, etc. containing
 * Vectorial data-structures.
 *
 * A common example would be \c std::vector< Math::Vector3d >
 * that is altered via \c std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __UBITRACK_MATH_BLAS_LEVEL_1_H__
#define __UBITRACK_MATH_BLAS_LEVEL_1_H__

#include "Util/vector_traits.h"
#include <utUtil/StaticAssert.h>


#include <cmath> // std::fabs, std::sqrt
#include <assert.h>

namespace Ubitrack { namespace Math {


/**
 * @ingroup math functor
 * @internal
 * Functor class to calculate the inner product 
 * of two vectors using a recursive implementation.
 */
struct InnerProduct	
{
public:
	
	/**
	 * @ingroup math
	 * Calculates the inner product of two vectors.
	 *
	 * @tparam VecType type of vector
	 * @param vec1 the 1st input vector 
	 * @param vec2 the 2nd input vector 
	 * @return inner product of the two vectors
	 */
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec1, const VecType& vec2 ) const
	{
		return this->operator()( vec1, vec2, typename Math::Util::vector_traits< VecType >::storage_category() );
	}
	
	/**
	 * @ingroup math
	 * Calculates the inner product of two vectors.
	 *
	 * Overloaded function for a single vector.
	 *
	 * @tparam VecType type of vector
	 * @param vec the input vector 
	 * @return inner product of the vector (== squared distance)
	 */
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec ) const
	{
		return this->operator()( vec, vec );
	}

protected:

	/// @internal function for fixed storage vectors, uses the internal functor
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec1, const VecType& vec2, const Math::Util::fixed_storage_tag ) const
	{
		UBITRACK_STATIC_ASSERT( ( Math::Util::has_fixed_storage< VecType >::value ), VECTOR_TYPE_OF_FIXED_STORAGE_CATEGORY_EXPECTED );
		typedef typename Math::Util::vector_traits< VecType >::value_type value_type;
		
		return inner_product_impl< value_type, Math::Util::vector_traits< VecType >::size >()( vec1, vec2, 0 );
	}
	
	/// @internal function for dynamic storage vectors, calculates the result using a loop
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec1, const VecType& vec2, const Math::Util::dynamic_storage_tag ) const
	{
	
		UBITRACK_STATIC_ASSERT( ( Math::Util::has_dynamic_storage< VecType >::value ), VECTOR_TYPE_OF_DYNAMIC_STORAGE_CATEGORY_EXPECTED );
		typedef typename Math::Util::vector_traits< VecType >::size_type size_type;
		typedef typename Math::Util::vector_traits< VecType >::value_type value_type;	
		
		const size_type size1 = Math::Util::vector_traits< VecType >::size( vec1 );
		const size_type size2 = Math::Util::vector_traits< VecType >::size( vec2 );
		assert( size1 == size2 );
		
		value_type result( 0 );
		for( size_type i = 0; i<size1; ++i )
			result += ( vec1[ i ] * vec2[ i ] );
		
		return result;
	}
	/**
	 * @ingroup math
	 * Internal functor that implements the calculation of the
	 * inner product of two vectors, also known as dot product.
	 *
	 * The implementation uses recursive functional programming,
	 * that optimizes calculations at compile time.
	 *
	 * @tparam T builtin-type of vector ( e.g. \c double or \c float )
	 * @tparam N dimension of vector (typically 2 or 3)
	 */
	template< typename T, std::size_t N >
	struct inner_product_impl
	{
		template< typename VecType >
		T operator() ( const VecType& vec1, const VecType& vec2, const T sqsum ) const
		{
			const T sq = sqsum + ( vec1[ N-1 ] * vec2[ N-1 ] );
			return inner_product_impl< T, N-1 >()( vec1, vec2, sq );
		}
	};

	/**
	 * Partial specialization of functor for final return value.
	 *
	 * @tparam T builtin-type of vector ( e.g. \c double or \c float )
	 */
	template< typename T >
	struct inner_product_impl< T, 0 >
	{
		template< typename VecType >
		T operator() ( const VecType&, const VecType& , const T sqsum ) const
		{
			return sqsum;
		}
	};
};

/**
 * @ingroup math functor
 * @internal
 * Functor class to calculate the 1-norm
 * (aka Manhatten norm) of a vector using a recursive
 * implementation.
 */
struct Norm_1
{
public:
	/**
	 * @ingroup math
	 * Calculates the 1-norm of a single vector.
	 *
	 * @tparam VecType type of vector
	 * @param vec the input vector 
	 * @return 1-norm of the vector
	 */
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec ) const
	{
		return this->operator()( vec, typename Math::Util::vector_traits< VecType >::storage_category() );
	}
	
protected:

	/// @internal function for fixed storage vectors, uses the internal functor
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec, const Math::Util::fixed_storage_tag ) const
	{
		UBITRACK_STATIC_ASSERT( ( Math::Util::has_fixed_storage< VecType >::value ), VECTOR_TYPE_OF_FIXED_STORAGE_CATEGORY_EXPECTED );
		typedef typename Math::Util::vector_traits< VecType >::value_type value_type;
		
		return norm_1_impl< value_type, Math::Util::vector_traits< VecType >::size >()( vec, 0 );
	}
	
	/// @internal function for dynamic storage vectors, calculates the result using a loop
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec, const Math::Util::dynamic_storage_tag ) const
	{
	
		UBITRACK_STATIC_ASSERT( ( Math::Util::has_dynamic_storage< VecType >::value ), VECTOR_TYPE_OF_DYNAMIC_STORAGE_CATEGORY_EXPECTED );
		typedef typename Math::Util::vector_traits< VecType >::size_type size_type;
		typedef typename Math::Util::vector_traits< VecType >::value_type value_type;	
		
		const size_type size = Math::Util::vector_traits< VecType >::size( vec );
		
		value_type result( 0 );
		for( size_type i = 0; i<size; ++i )
			result += std::fabs( vec[ i ] );
		
		return result;
	}
	
	/**
	 * @ingroup math
	 * Internal functor that implements the calculation of the
	 * 1-norm of a single vectors, also known as Manhattan norm.
	 *
	 * The implementation uses recursive functional programming,
	 * that optimizes calculations at compile time.
	 *
	 * @tparam T builtin-type of vector ( e.g. \c double or \c float )
	 * @tparam N dimension of vector (typically 2 or 3)
	 */
	template< typename T, std::size_t N >
	struct norm_1_impl
	{
		template< typename VecType >
		T operator() ( const VecType& vec, const T sqsum ) const
		{
			const T sq = sqsum + std::fabs( vec[ N-1 ] );
			return norm_1_impl< T, N-1 >()( vec, sq );
		}
	};

	/**
	 * Partial specialization of functor for final return value.
	 *
	 * @tparam T builtin-type of vector ( e.g. \c double or \c float )
	 */
	template< typename T >
	struct norm_1_impl< T, 0 >
	{
		template< typename VecType >
		T operator() ( const VecType&, const T sqsum ) const
		{
			return sqsum;
		}
	};
};

/**
 * @ingroup math functor
 * @internal
 * Functor class to calculate the 2-norm of a vector
 * (aka Euclidean norm) using a recursive
 * implementation.
 */
struct Norm_2
{
public:
	/**
	 * @ingroup math
	 * Calculates the 2-norm of a single vector.
	 *
	 * @tparam VecType type of vector
	 * @param vec the input vector 
	 * @return 2-norm of the vector
	 */
	template< typename VecType >
	typename Math::Util::vector_traits< VecType >::value_type operator() ( const VecType& vec ) const
	{
		return std::sqrt( InnerProduct().operator()< VecType >( vec, vec ) );
	}
};

/**
 * @ingroup math
 * @brief A function that calculates the inner product of two vectors.
 * 
 * This function calculates the inner product of two vectors \f$ u \f$
 * and \f$ v \f$ , each one consisting of n elements, as
 * \f$ u^T \cdot v = \sum_{i=1}^{n} u_i \cdot v_i \f$
 * , which is also known as the dot product or scalar product.
 *
 * This function template wraps a call to the InnerProduct
 * functor. The input vectors for calculating the inner product can be
 * of any dimension. 
 *
 * @tparam VecType the type of the input vectors. The type must be the same for both the input vectors.
 * @param vec1 the \b 1st input vector
 * @param vec2 the \b 2nd input vector
 * @return the innner product of the two vectors
 */
template< typename VecType >
inline typename Math::Util::vector_traits< VecType >::value_type inner_product( const VecType& vec1, const VecType& vec2 )
{
	return InnerProduct()( vec1, vec2 );
}

/**
 * @ingroup math
 * @brief A function that calculates the 1-norm of a vector.
 * 
 * This function calculates the 1-norm of a vector \f$ v \f$ with
 * n elements as \f$ ||v||_1 = \sum_{i=1}^{n} |v_i| \f$
 * , which is also known as the Manhattan norm of a vector.
 *
 * This function template wraps a call to the Norm_1 functor.
 * The input vector for calculating the 1-norm can be
 * of any dimension.
 *
 * @tparam VecType the type of the input vector.
 * @param vec the input vector
 * @return the 1-norm of the input vector
 */
template< typename VecType >
inline typename Math::Util::vector_traits< VecType >::value_type norm_1( const VecType& vec )
{
	return Norm_1().operator()< VecType >( vec );
}

/**
 * @ingroup math
 * @brief A function that calculates the 2-norm of a vector.
 * 
 * This function calculates the 2-norm of a vector \f$ v \f$ with
 * n elements as \f$ ||v||_2 = \sqrt{ \sum_{i=1}^{n} v_i^2 } \f$
 * , which is also known as the Euclidean norm of a vector.
 *
 * This function template wraps a call to the Norm_2 functor.
 * The input vector for calculating the 2-norm can be
 * of any dimension.
 *
 * @tparam VecType the type of the input vector.
 * @param vec the input vector
 * @return the 2-norm of the input vector
 */
template< typename VecType >
inline typename Math::Util::vector_traits< VecType >::value_type norm_2( const VecType& vec )
{
	return Norm_2().operator()< VecType >( vec );
}

} } // namespace Ubitrack::Math

#endif //__UBITRACK_MATH_BLAS_LEVEL_1_H__

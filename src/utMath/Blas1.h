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

 
 
#ifndef __H__BLAS_LEVEL_1__
#define __H__BLAS_LEVEL_1__

#include "vector_traits.h"

namespace Ubitrack { namespace Math {


/**
 * @ingroup math
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
	typename Math::vector_traits< VecType >::value_type operator() ( const VecType& vec1, const VecType& vec2 ) const
	{
		UBITRACK_STATIC_ASSERT( ( Math::has_fixed_storage< VecType >::value ), NEED_VECTOR_OF_FIXED_STORAGE_TYPE );
		typedef Math::vector_traits< VecType >::value_type value_type;
		
		return inner_product_impl< value_type, Math::vector_traits< VecType >::size >()( vec1, vec2, 0 );
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
	typename Math::vector_traits< VecType >::value_type operator() ( const VecType& vec ) const
	{
		return this->operator()( vec, vec );
	}

protected:
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
 * @ingroup math
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
	typename Math::vector_traits< VecType >::value_type operator() ( const VecType& vec ) const
	{
		UBITRACK_STATIC_ASSERT( ( Math::has_fixed_storage< VecType >::value ), NEED_VECTOR_OF_FIXED_STORAGE_TYPE );
		typedef Math::vector_traits< VecType >::value_type value_type;
		
		return norm_1_impl< value_type, Math::vector_traits< VecType >::size >()( vec, 0 );
	}
	
protected:
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
 * @ingroup math
 * Functor class to calculate the Euclidean norm
 * (aka 2-norm) of a vector using a recursive
 * implementation.
 */
struct Norm_2
{
protected:
	static const InnerProduct dotterer;
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
	typename Math::vector_traits< VecType >::value_type operator() ( const VecType& vec ) const
	{
		UBITRACK_STATIC_ASSERT( ( Math::has_fixed_storage< VecType >::value ), NEED_VECTOR_OF_FIXED_STORAGE_TYPE );
		typedef Math::vector_traits< VecType >::value_type value_type;
		
		return std::sqrt( dotterer.operator()( vec ) );
	}
};

} } // namespace Ubitrack::Math

#endif //__H__BLAS_LEVEL_1__

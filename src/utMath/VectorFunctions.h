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
 * Functors and functions for common higher level vector
 * operations. The functions inside this file represent
 * common functions, other than BLAS level 1 operations.
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


#ifndef __UBTRACK_MATH_VECTOR_FUNCTIONS_H__
#define __UBTRACK_MATH_VECTOR_FUNCTIONS_H__

// Ubitrack
#include "Blas1.h"

// std
#include <iterator> // iterator_traits

namespace Ubitrack { namespace Math {

/**
 * @ingroup math functor
 * @internal
 * Functor class to normalize a vector by its Euclidean norm.
 *
 * @tparam VecType type of vector (e.g. \c Math::Vector3d or \c Math::Vector3f )
 */
template< typename VecType >
struct Normalize
{
public:
	/**
	 * @ingroup math
	 * @internal
	 * Normalizes a Math::Vector such that his 
	 * length equals to one.
	 *
	 *
	 * @param vec the input vector 
	 * @return normalized vector
	 */
	VecType operator() ( const VecType& vec ) const
	{
		const typename Math::Util::vector_traits< VecType >::value_type norm ( Norm_2()( vec ) );
		return vec / norm;
	}
};

/**
 * @ingroup math functor
 * @internal
 * Functor class to calculate the Euclidean distance 
 * between two vectors.
 *
 * @tparam VecType type of vector (e.g. \c Math::Vector3d or \c Math::Vector3f )
 */
template< typename VecType >
struct Distance
{
public:

	/**
	 * @ingroup math
	 * @internal
	 * Calculates the Euclidean distance of a vector to the origin.
	 *
	 * @param vec the input vector 
	 * @return length of the vector (Euclidean distance between the origin and the vector)
	 */
	typename Math::Util::vector_traits< VecType >::value_type operator()( const VecType& vec ) const
	{
		return Ubitrack::Math::Norm_2()( vec );
	}

	/**
	 * @ingroup math
	 * @internal
	 * Calculates the Euclidean distance between two vectors.
	 *
	 * @param vec1 the 1st input vector 
	 * @param vec2 the 2nd input vector 
	 * @return Euclidean distance between the two vectors
	 */
	typename Math::Util::vector_traits< VecType >::value_type operator()( const VecType& vec1, const VecType& vec2 ) const
	{
		const VecType vec ( vec1 - vec2 );
		return this->operator()( vec );
	}
};


/**
 * @ingroup math functor
 * @internal
 * Functor class to calculate the squared Euclidean distance 
 * between two vectors.
 *
 * @tparam VecType type of vector (e.g. \c Math::Vector3d or \c Math::Vector3f )
 */
template< typename VecType >
struct SquaredDistance
{
public:

	/**
	 * @ingroup math
	 * @internal
	 * Calculates the squared Euclidean distance of a vector to the origin.
	 *
	 * @param vec the input vector 
	 * @return squared length of the vector (squared Euclidean distance between the origin and the vector)
	 */
	typename Math::Util::vector_traits< VecType >::value_type operator()( const VecType& vec ) const
	{
		return Ubitrack::Math::InnerProduct()( vec );
	}

	/**
	 * @ingroup math
	 * @internal
	 * Calculates the squared Euclidean distance between two vectors.
	 *
	 * @param vec1 the 1st input vector 
	 * @param vec2 the 2nd input vector 
	 * @return squared Euclidean distance between the two vectors
	 */
	typename Math::Util::vector_traits< VecType >::value_type operator()( const VecType& vec1, const VecType& vec2 ) const
	{
		const VecType vec ( vec1 - vec2 );
		return this->operator()( vec );
	}
};

/**
 * @ingroup math
 * @brief A function that calculates the Euclidean distance between two vectors.
 * 
 * This function calculates the Euclidean distance \f$ d \f$ between two
 * vectors \f$ v_1 \f$ and \f$ v_2 \f$ as \f$  d = ||v_1 - v_2||_2 \f$ .
 *
 * This function template wraps a call to the \c Distance functor.
 * The input vector for calculating the Euclidean distance can be of any
 * dimension.
 *
 * @tparam VecType the type of the input vector.
 * @param vec1 the 1st input vector
 * @param vec2 the 2nd input vector
 * @return the Euclidean distance between the two vectors.
 */
template< typename VecType >
inline typename Math::Util::vector_traits< VecType >::value_type distance( const VecType& vec1, const VecType& vec2 )
{
	return Distance< VecType >().operator()( vec1, vec2 );
}

/**
 * @ingroup math
 * @brief A function that calculates the Euclidean distance between vectors and the origin.
 * 
 * This function calculates the Euclidean distance \f$ d_i \f$ between 
 * vectors \f$ v_i \f$ and the origin \f$  d_i = ||v_i ||_2 \f$ .
 *
 * This function template wraps a call to the \c std::transform with
 * \c Distance functor as the unary function.
 * The input vectors for calculating the Euclidean Distance can be
 * of any dimension.
 *
 * @tparam InputIterator type of forward iterator to container of input vectors
 * @tparam OutputIterator type of forward iterator to container of Euclidean distances
 * @param iBegin \c iterator pointing to first element in the input container/storage class of the \b vectors ( usually \c begin() )
 * @param iEnd \c iterator pointing behind the last element in the input container/storage class of the \b vectors ( usually \c end() )
 * @param iOut output \c iterator pointing to first element in container/storage class for storing the Euclidean distances as a result ( usually \c begin() or \c std::back_inserter(container) )
 */
template< typename InputIterator, typename OutputIterator >
inline void distance( const InputIterator iBegin, const InputIterator iEnd, OutputIterator iOut )
{
	std::transform( iBegin, iEnd, iOut, Distance< typename std::iterator_traits< InputIterator >::value_type >() );
}

/**
 * @ingroup math
 * @brief A function that calculates the Euclidean distance between two vectors.
 * 
 * This function calculates the Euclidean distance \f$ d_i \f$ between 
 * two vectors \f$ v1_i \f$ and \f$ v2_i \f$ as \f$  d_i = ||v1_i - v2_i ||_2 \f$ .
 *
 * This function template wraps a call to the \c std::transform with
 * \c Distance functor as the binary function.
 * The input vectors for calculating the Euclidean Distance can be
 * of any dimension.
 *
 * @tparam InputIterator type of forward iterator to container of input vectors
 * @tparam OutputIterator type of forward iterator to container of Euclidean distances
 * @param iBegin1 \c iterator pointing to first element in the 1st input container/storage class of the \b vectors ( usually \c begin() )
 * @param iEnd1 \c iterator pointing behind the last element in the 1st input container/storage class of the \b vectors ( usually \c end() )
 * @param iBegin2 \c iterator pointing to first element in the 2nd input container/storage class of the \b vectors ( usually \c begin() )
 * @param iOut output \c iterator pointing to first element in container/storage class for storing the Euclidean distances as a result ( usually \c begin() or \c std::back_inserter(container) )
 */
template< typename InputIterator, typename OutputIterator >
inline void distance( const InputIterator iBegin1, const InputIterator iEnd1, const InputIterator iBegin2, OutputIterator iOut )
{
	std::transform( iBegin1, iEnd1, iBegin2, iOut, Distance< typename std::iterator_traits< InputIterator >::value_type >() );
}

/**
 * @ingroup math
 * @brief A function that calculates the squared Euclidean distance between two vectors.
 * 
 * This function calculates the squared Euclidean distance \f$ d \f$ between two
 * vectors \f$ v_1 \f$ and \f$ v_2 \f$ as \f$  d = ||v_1 - v_2||_2^2 \f$ .
 *
 * This function template wraps a call to the \c SquaredDistance functor.
 * The input vector for calculating the squared Euclidean distance can be of any
 * dimension.
 *
 * @tparam VecType the type of the input vector.
 * @param vec1 the 1st input vector
 * @param vec2 the 2nd input vector
 * @return the squared Euclidean distance between the two vectors.
 */
template< typename VecType >
inline typename Math::Util::vector_traits< VecType >::value_type squared_distance( const VecType& vec1, const VecType& vec2 )
{
	return SquaredDistance< VecType >().operator()( vec1, vec2 );
}

/**
 * @ingroup math
 * @brief A function that calculates the squared Euclidean distance between vectors and the origin.
 * 
 * This function calculates the Euclidean distance \f$ d_i \f$ between 
 * vectors \f$ v_i \f$ and the origin \f$  d_i = || v_i ||_2^2 \f$ .
 *
 * This function template wraps a call to the \c std::transform with
 * \c SquaredDistance functor as the unary function.
 * The input vectors for calculating the squared Euclidean Distance can be
 * of any dimension.
 *
 * @tparam InputIterator type of forward iterator to container of input vectors
 * @tparam OutputIterator type of forward iterator to container of squared Euclidean distances
 * @param iBegin \c iterator pointing to first element in the input container/storage class of the \b vectors ( usually \c begin() )
 * @param iEnd \c iterator pointing behind the last element in the input container/storage class of the \b vectors ( usually \c end() )
 * @param iOut output \c iterator pointing to first element in container/storage class for storing the squared Euclidean distances as a result ( usually \c begin() or \c std::back_inserter(container) )
 */
template< typename InputIterator, typename OutputIterator >
inline void squared_distance( const InputIterator iBegin, const InputIterator iEnd, OutputIterator iOut )
{
	std::transform( iBegin, iEnd, iOut, SquaredDistance< typename std::iterator_traits< InputIterator >::value_type >() );
}

/**
 * @ingroup math
 * @brief A function that calculates the squared Euclidean distance between two vectors.
 * 
 * This function calculates the Euclidean distance \f$ d_i \f$ between 
 * two vectors \f$ v1_i \f$ and \f$ v2_i \f$ as \f$  d_i = || v1_i - v2_i ||_2^2 \f$ .
 *
 * This function template wraps a call to the \c std::transform with
 * \c SquaredDistance functor as the binary function.
 * The input vectors for calculating the squared Euclidean distance can be
 * of any dimension.
 *
 * @tparam InputIterator type of forward iterator to container of input vectors
 * @tparam OutputIterator type of forward iterator to container of squared Euclidean distances
 * @param iBegin1 \c iterator pointing to first element in the 1st input container/storage class of the \b vectors ( usually \c begin() )
 * @param iEnd1 \c iterator pointing behind the last element in the 1st input container/storage class of the \b vectors ( usually \c end() )
 * @param iBegin2 \c iterator pointing to first element in the 2nd input container/storage class of the \b vectors ( usually \c begin() )
 * @param iOut output \c iterator pointing to first element in container/storage class for storing the squared Euclidean distances as a result ( usually \c begin() or \c std::back_inserter(container) )
 */
template< typename InputIterator, typename OutputIterator >
inline void squared_distance( const InputIterator iBegin1, const InputIterator iEnd1, const InputIterator iBegin2, OutputIterator iOut )
{
	std::transform( iBegin1, iEnd1, iBegin2, iOut, SquaredDistance< typename std::iterator_traits< InputIterator >::value_type >() );
}

/**
 * @ingroup math
 * @brief A function that normalizes a vector by it's 2-norm.
 * 
 * This function normalizes a vector \f$ v \f$ with n elements
 * such that each single element is divided by the 2-norm of
 * the vector as \f$  v = v / ||v||_2 \f$ .
 *
 * This function template wraps a call to the NormalizeVector functor.
 * The input vector for calculating the normalized vector can be
 * of any dimension.
 *
 * @tparam VecType the type of the input vector.
 * @param vec the input vector
 * @return the normalized input vector
 */
template< typename VecType >
inline VecType normalize( const VecType& vec )
{
	return Normalize< VecType >().operator()( vec );
}

/**
 * @ingroup math
 * @brief A function that normalizes vectors by their 2-norm.
 * 
 * This function normalizes vectors \f$ v_i \f$ with n elements
 * such that each single element of a vector \f$ v_i \f$
 *  is divided by its 2-norm of as \f$  v_i = v_i / ||v_i||_2 \f$ .
 *
 * This function template wraps a call to the \c std::transform with
 * \c NormalizeVector functor as the unary function.
 * The input vectors for calculating the normalized vectors can be
 * of any dimension.
 *
 * @tparam InputIterator type of forward iterator to container of input and output vectors
 * @param iBegin \c iterator pointing to first element in the input container/storage class of the \b vectors ( usually \c begin() )
 * @param iEnd \c iterator pointing behind the last element in the input container/storage class of the \b vectors ( usually \c end() )
 * @param iOut output \c iterator pointing to first element in container/storage class for storing the transformed vectors as a result ( usually \c begin() or \c std::back_inserter(container) )
 */
template< typename ForwardIterator >
inline void normalize( const ForwardIterator iBegin, const ForwardIterator iEnd, ForwardIterator iOut )
{
	std::transform( iBegin, iEnd, iOut, Normalize< typename std::iterator_traits< ForwardIterator >::value_type >() );
}

} } // namespace Ubitrack::Math

#endif //__UBTRACK_MATH_VECTOR_FUNCTIONS_H__

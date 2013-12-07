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
 * std::vector of Math::Vectors< T, N > using std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

 
 
#ifndef __H__VECTOR_FUNCTORS__
#define __H__VECTOR_FUNCTORS__


#include "Vector2Functors.h"
#include "Vector3Functors.h"
#include "VectorNFunctors.h"

#include <boost/numeric/ublas/matrix.hpp>

#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>



namespace Ubitrack { namespace Math { namespace Functors {


/**
 * @ingroup math
 * Transforms the Math::Vector< T, 3 > with a 3x4 transformation matrix.
 */
template< typename T >
struct TransformVector
{
protected: 
	const Ubitrack::Math::Matrix< T, 3, 4 > m_transformation;
	
public:

	TransformVector( const Ubitrack::Math::Matrix< T, 3, 4 > &transformation )
		: m_transformation( transformation )
	{};
	
	TransformVector( const Ubitrack::Math::Pose &pose )
		: m_transformation( pose )
	{};
	
	TransformVector( const Ubitrack::Math::Quaternion &rotation, const Math::Vector< T, 3 > &translation )
		: m_transformation( rotation, translation )
	{};

	Ubitrack::Math::Vector< T, 3 > operator() ( const Ubitrack::Math::Vector< T, 2 > &vec ) const
	{
		const T e1 = m_transformation( 0, 0 ) * vec( 0 ) + m_transformation( 0, 1 ) * vec( 1 ) + m_transformation( 0, 2 );
		const T e2 = m_transformation( 1, 0 ) * vec( 0 ) + m_transformation( 1, 1 ) * vec( 1 ) + m_transformation( 1, 2 );
		const T e3 = m_transformation( 2, 0 ) * vec( 0 ) + m_transformation( 2, 1 ) * vec( 1 ) + m_transformation( 2, 2 );
		return Ubitrack::Math::Vector< T, 3 > ( e1, e2, e3 );
	}
	
	Ubitrack::Math::Vector< T, 3 > operator() ( const Ubitrack::Math::Vector< T, 3 > &vec ) const
	{
		const T e1 = m_transformation( 0, 0 ) * vec( 0 ) + m_transformation( 0, 1 ) * vec( 1 ) + m_transformation( 0, 2 ) * vec( 2 ) + m_transformation( 0, 3 );
		const T e2 = m_transformation( 1, 0 ) * vec( 0 ) + m_transformation( 1, 1 ) * vec( 1 ) + m_transformation( 1, 2 ) * vec( 2 ) + m_transformation( 1, 3 );
		const T e3 = m_transformation( 2, 0 ) * vec( 0 ) + m_transformation( 2, 1 ) * vec( 1 ) + m_transformation( 2, 2 ) * vec( 2 ) + m_transformation( 2, 3 );
		return Ubitrack::Math::Vector< T, 3 > ( e1, e2, e3 );
	}
	
	Ubitrack::Math::Vector< T, 3 > operator() ( const Math::Vector< T, 4 > &vec ) const
	{
		const T e1 = m_transformation( 0, 0 ) * vec( 0 ) + m_transformation( 0, 1 ) * vec( 1 ) + m_transformation( 0, 2 ) * vec( 2 ) + m_transformation( 0, 3 ) * vec( 3 );
		const T e2 = m_transformation( 1, 0 ) * vec( 0 ) + m_transformation( 1, 1 ) * vec( 1 ) + m_transformation( 1, 2 ) * vec( 2 ) + m_transformation( 1, 3 ) * vec( 3 );
		const T e3 = m_transformation( 2, 0 ) * vec( 0 ) + m_transformation( 2, 1 ) * vec( 1 ) + m_transformation( 2, 2 ) * vec( 2 ) + m_transformation( 2, 3 ) * vec( 3 );
		return Ubitrack::Math::Vector< T, 3 > ( e1, e2, e3 );
	}
};


/**
 * @ingroup math
 * Projects a Math::Vector by a given 3x4 projection matrix.
 * several overloads for 2, 3 and 4 vector are provided.
 */
template< typename T >
struct ProjectVector
{
protected: 
	const Ubitrack::Math::Matrix< T, 3, 4 > m_projection;
	
public:

	ProjectVector( const Ubitrack::Math::Matrix< T, 3, 4 > &projection )
		: m_projection( projection )
	{};
	
	ProjectVector( const Ubitrack::Math::Matrix< T, 3, 3 > &projection, const Ubitrack::Math::Pose &pose )
		: m_projection( boost::numeric::ublas::prod( projection, Ubitrack::Math::Matrix< T, 3, 4 >( pose ) ) )
	{};
	
	ProjectVector( const Ubitrack::Math::Matrix< T, 3, 3 > &projection, const Ubitrack::Math::Quaternion &rotation, const Ubitrack::Math::Vector< T, 3 > &translation )
		: m_projection( boost::numeric::ublas::prod( projection, Ubitrack::Math::Matrix< T, 3, 4 >( rotation, translation ) ) )
	{};

	Ubitrack::Math::Vector< T, 2 > operator() ( const Ubitrack::Math::Vector< T, 2 > &vec ) const
	{
		const T e1 = m_projection( 0, 0 ) * vec( 0 ) + m_projection( 0, 1 ) * vec( 1 ) + m_projection( 0, 3 );
		const T e2 = m_projection( 1, 0 ) * vec( 0 ) + m_projection( 1, 1 ) * vec( 1 ) + m_projection( 1, 3 );
		const T e3 = m_projection( 2, 0 ) * vec( 0 ) + m_projection( 2, 1 ) * vec( 1 ) + m_projection( 2, 3 );
		return Ubitrack::Math::Vector< T, 2 > ( e1/e3, e2/e3 );
	}
	
	Ubitrack::Math::Vector< T, 2 > operator() ( const Ubitrack::Math::Vector< T, 3 > &vec ) const
	{
		const T e1 = m_projection( 0, 0 ) * vec( 0 ) + m_projection( 0, 1 ) * vec( 1 ) + m_projection( 0, 2 ) * vec( 2 ) + m_projection( 0, 3 );
		const T e2 = m_projection( 1, 0 ) * vec( 0 ) + m_projection( 1, 1 ) * vec( 1 ) + m_projection( 1, 2 ) * vec( 2 ) + m_projection( 1, 3 );
		const T e3 = m_projection( 2, 0 ) * vec( 0 ) + m_projection( 2, 1 ) * vec( 1 ) + m_projection( 2, 2 ) * vec( 2 ) + m_projection( 2, 3 );
		return Ubitrack::Math::Vector< T, 2 > ( e1/e3, e2/e3 );
	}
	
	Ubitrack::Math::Vector< T, 2 > operator() ( const Ubitrack::Math::Vector< T, 4 > &vec ) const
	{
		const T e1 = m_projection( 0, 0 ) * vec( 0 ) + m_projection( 0, 1 ) * vec( 1 ) + m_projection( 0, 2 ) * vec( 2 ) + m_projection( 0, 3 ) * vec( 3 );
		const T e2 = m_projection( 1, 0 ) * vec( 0 ) + m_projection( 1, 1 ) * vec( 1 ) + m_projection( 1, 2 ) * vec( 2 ) + m_projection( 1, 3 ) * vec( 3 );
		const T e3 = m_projection( 2, 0 ) * vec( 0 ) + m_projection( 2, 1 ) * vec( 1 ) + m_projection( 2, 2 ) * vec( 2 ) + m_projection( 2, 3 ) * vec( 3 );
		return Math::Vector< T, 2 > ( e1/e3, e2/e3 );
	}
};


/**
 * @ingroup math
 * Functor class to calculate the 1st norm
 * of a vector.
 *
 * Recursive implementation that can be applied to 
 * a Math::Vector of arbitrary dimension N.
 *
 * @tparam N dimension of vector
 * @tparam T builtin-type of vector
 */
template< typename T, std::size_t N >
struct Norm_1
{
public:
	/**
	 * @ingroup math
	 * Calculates the 1st norm of a Math::Vector.
	 *
	 * @param vec the vector 
	 * @return norm of the N-dimensional vector
	 */
	T operator() ( const Math::Vector< T, N >& vec ) const
	{
		return norm1_impl< N, Math::Vector< T, N > >()( vec, 0 );
	}

protected:
	/**
	 * @ingroup math
	 * Internal functor for recursive implementation.
	 *
	 */
	template< std::size_t I, typename VecType >
	struct norm1_impl
	{
		T operator() ( const VecType& vec, const T sqsum ) const
		{
			const T sq = sqsum + ( vec[ I-1 ] * vec[ I-1 ] );
			return norm1_impl< I-1, VecType >()( vec, sq );
		}
	};

	/**
	 * @ingroup math
	 * Partial specialization of functor for final element.
	 *
	 */
	template< typename VecType >
	struct norm1_impl< 0, VecType >
	{
		T operator() ( const VecType& vec, const T sqsum ) const
		{
			return sqsum;
		}
	};
};


/**
 * @ingroup math
 * Functor Class to calcualte the Euclidean distance
 * of a vector.
 *
 * Uses internally the functor that calcualtes the
 * 1st norm of a vector and can be applied to
 * a Math::Vector of arbitrary dimension N.
 *
 * @tparam N dimension of vector
 * @tparam T builtin-type of vector
 */
template< typename T, std::size_t N >
struct Norm_2
{
public:
	/**
	 * @ingroup math
	 * Calculates the length of a Math::Vector.
	 *
	 * @param vec the vector 
	 * @return norm of the N-dimensional vector
	 */
	T operator() ( const Math::Vector< T, N >& vec ) const
	{
		return std::sqrt( Norm_1< T, N >().operator()( vec ) );
	}
};


/**
 * @ingroup math
 * Functor Class to the calculate a normalized vector.
 *
 * Recursive implementation for Vectors of
 * arbitrary dimension N.
 *
 * @tparam N dimension of vector
 * @tparam T builtin-type of vector
 */
template< typename T, std::size_t N >
struct NormalizeVector
	: public std::unary_function< Math::Vector< T, N >, Math::Vector< T, N > >
{
public:
	/**
	 * @ingroup math
	 * Normalizes a Math::Vector such that his 
	 * length equals to one.
	 *
	 * @param vec the N-dimensional vector 
	 * @return norm of the N-dimensional vector
	 */
	Math::Vector< T, N > operator() ( const Math::Vector< T, N >& vec ) const
	{
		const T norm ( Norm_2< T, N >()( vec ) );
		return vec * ( 1./ norm );
	}
};


/**
 * @ingroup math
 * Functor Class to the calculate the inner product 
 * of two vectors.
 *
 * Recursive implementation for Vectors of
 * arbitrary dimension N.
 *
 * @tparam N dimension of vector
 * @tparam T builtin-type of vector
 */
template< std::size_t N, typename T >
struct InnerProduct
	: public std::binary_function< Math::Vector< T, N >, Math::Vector< T, N >, T >
{
public:
	/**
	 * @ingroup math
	 * Calculates the inner product of two vectors.
	 *
	 * @param vec1 the 1st vector 
	 * @param vec2 the 2nd vector 
	 * @return inner product of the two vectors
	 */
	T operator() ( const Math::Vector< T, N >& vec1, const Math::Vector< T, N >& vec2 ) const
	{
		return inner_product_impl< N, Math::Vector< T, N > >()( vec1, vec2, 0 );
	}

protected:
	/**
	 * @ingroup math
	 * Internal functor for recursive implementation 
	 * of dot product.
	 *
	 */
	template< std::size_t I, typename VecType >
	struct inner_product_impl
	{
		T operator() ( const VecType& vec1, const VecType& vec2, const T sqsum ) const
		{
			const T sq = sqsum + ( vec1[ I-1 ] * vec2[ I-1 ] );
			return inner_product_impl< I-1, VecType >()( vec1, vec2, sq );
		}
	};

	/**
	 * @ingroup math
	 * Partial specialization of functor for final element.
	 *
	 */
	template< typename VecType >
	struct inner_product_impl< 0, VecType >
	{
		T operator() ( const VecType&, const VecType& , const T sqsum ) const
		{
			return sqsum;
		}
	};
};


} } } // namespace Ubitrack::Math::Functors

#endif //__H__VECTOR_FUNCTORS__

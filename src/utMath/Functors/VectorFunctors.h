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
 * Transforms the Math::Vector< 3, T > with a 3x4 transformation matrix.
 */
template< typename T >
struct TransformVector
{
protected: 
	const Ubitrack::Math::Matrix< 3, 4, T > m_transformation;
	
public:

	TransformVector( const Ubitrack::Math::Matrix< 3, 4, T > &transformation )
		: m_transformation( transformation )
	{};
	
	TransformVector( const Ubitrack::Math::Pose &pose )
		: m_transformation( pose )
	{};
	
	TransformVector( const Ubitrack::Math::Quaternion &rotation, const Math::Vector< 3, T > &translation )
		: m_transformation( rotation, translation )
	{};

	Ubitrack::Math::Vector< 3, T > operator() ( const Ubitrack::Math::Vector< 2, T > &vec ) const
	{
		const T e1 = m_transformation( 0, 0 ) * vec( 0 ) + m_transformation( 0, 1 ) * vec( 1 ) + m_transformation( 0, 2 );
		const T e2 = m_transformation( 1, 0 ) * vec( 0 ) + m_transformation( 1, 1 ) * vec( 1 ) + m_transformation( 1, 2 );
		const T e3 = m_transformation( 2, 0 ) * vec( 0 ) + m_transformation( 2, 1 ) * vec( 1 ) + m_transformation( 2, 2 );
		return Ubitrack::Math::Vector< 3, T > ( e1, e2, e3 );
	}
	
	Ubitrack::Math::Vector< 3, T > operator() ( const Ubitrack::Math::Vector< 3, T > &vec ) const
	{
		const T e1 = m_transformation( 0, 0 ) * vec( 0 ) + m_transformation( 0, 1 ) * vec( 1 ) + m_transformation( 0, 2 ) * vec( 2 ) + m_transformation( 0, 3 );
		const T e2 = m_transformation( 1, 0 ) * vec( 0 ) + m_transformation( 1, 1 ) * vec( 1 ) + m_transformation( 1, 2 ) * vec( 2 ) + m_transformation( 1, 3 );
		const T e3 = m_transformation( 2, 0 ) * vec( 0 ) + m_transformation( 2, 1 ) * vec( 1 ) + m_transformation( 2, 2 ) * vec( 2 ) + m_transformation( 2, 3 );
		return Ubitrack::Math::Vector< 3, T > ( e1, e2, e3 );
	}
	
	Ubitrack::Math::Vector< 3, T > operator() ( const Math::Vector< 4, T > &vec ) const
	{
		const T e1 = m_transformation( 0, 0 ) * vec( 0 ) + m_transformation( 0, 1 ) * vec( 1 ) + m_transformation( 0, 2 ) * vec( 2 ) + m_transformation( 0, 3 ) * vec( 3 );
		const T e2 = m_transformation( 1, 0 ) * vec( 0 ) + m_transformation( 1, 1 ) * vec( 1 ) + m_transformation( 1, 2 ) * vec( 2 ) + m_transformation( 1, 3 ) * vec( 3 );
		const T e3 = m_transformation( 2, 0 ) * vec( 0 ) + m_transformation( 2, 1 ) * vec( 1 ) + m_transformation( 2, 2 ) * vec( 2 ) + m_transformation( 2, 3 ) * vec( 3 );
		return Ubitrack::Math::Vector< 3, T > ( e1, e2, e3 );
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
	const Ubitrack::Math::Matrix< 3, 4, T > m_projection;
	
public:

	ProjectVector( const Ubitrack::Math::Matrix< 3, 4, T > &projection )
		: m_projection( projection )
	{};
	
	ProjectVector( const Ubitrack::Math::Matrix< 3, 3, T > &projection, const Ubitrack::Math::Pose &pose )
		: m_projection( boost::numeric::ublas::prod( projection, Ubitrack::Math::Matrix< 3, 4, T >( pose ) ) )
	{};
	
	ProjectVector( const Ubitrack::Math::Matrix< 3, 3, T > &projection, const Ubitrack::Math::Quaternion &rotation, const Ubitrack::Math::Vector< 3, T > &translation )
		: m_projection( boost::numeric::ublas::prod( projection, Ubitrack::Math::Matrix< 3, 4, T >( rotation, translation ) ) )
	{};

	Ubitrack::Math::Vector< 2, T > operator() ( const Ubitrack::Math::Vector< 2, T > &vec ) const
	{
		const T e1 = m_projection( 0, 0 ) * vec( 0 ) + m_projection( 0, 1 ) * vec( 1 ) + m_projection( 0, 2 );
		const T e2 = m_projection( 1, 0 ) * vec( 0 ) + m_projection( 1, 1 ) * vec( 1 ) + m_projection( 1, 2 );
		const T e3 = m_projection( 2, 0 ) * vec( 0 ) + m_projection( 2, 1 ) * vec( 1 ) + m_projection( 2, 2 );
		return Ubitrack::Math::Vector< 2, T > ( e1/e3, e2/e3 );
	}
	
	Ubitrack::Math::Vector< 2, T > operator() ( const Ubitrack::Math::Vector< 3, T > &vec ) const
	{
		const T e1 = m_projection( 0, 0 ) * vec( 0 ) + m_projection( 0, 1 ) * vec( 1 ) + m_projection( 0, 2 ) * vec( 2 ) + m_projection( 0, 3 );
		const T e2 = m_projection( 1, 0 ) * vec( 0 ) + m_projection( 1, 1 ) * vec( 1 ) + m_projection( 1, 2 ) * vec( 2 ) + m_projection( 1, 3 );
		const T e3 = m_projection( 2, 0 ) * vec( 0 ) + m_projection( 2, 1 ) * vec( 1 ) + m_projection( 2, 2 ) * vec( 2 ) + m_projection( 2, 3 );
		return Ubitrack::Math::Vector< 2, T > ( e1/e3, e2/e3 );
	}
	
	Ubitrack::Math::Vector< 2, T > operator() ( const Ubitrack::Math::Vector< 4, T > &vec ) const
	{
		const T e1 = m_projection( 0, 0 ) * vec( 0 ) + m_projection( 0, 1 ) * vec( 1 ) + m_projection( 0, 2 ) * vec( 2 ) + m_projection( 0, 3 ) * vec( 3 );
		const T e2 = m_projection( 1, 0 ) * vec( 0 ) + m_projection( 1, 1 ) * vec( 1 ) + m_projection( 1, 2 ) * vec( 2 ) + m_projection( 1, 3 ) * vec( 3 );
		const T e3 = m_projection( 2, 0 ) * vec( 0 ) + m_projection( 2, 1 ) * vec( 1 ) + m_projection( 2, 2 ) * vec( 2 ) + m_projection( 2, 3 ) * vec( 3 );
		return Math::Vector< 2, T > ( e1/e3, e2/e3 );
	}
};

/**
 * @ingroup math
 * Functor Class to the calculate the 1st norm
 * of a vector.
 *
 * Recursive implementation for a Math::Vector 
 * of arbitrary dimension N.
 */
struct Norm_1
{
public:
	/**
	 * @ingroup math
	 * Calculates the length of a Math::Vector.
	 *
	 * @tparam N dimension of vector
	 * @tparam T builtin-type of vector
	 * @param vec the vector 
	 * @return norm of the N-dimensional vector
	 */
	template< std::size_t N, typename T >
	T operator() ( const Math::Vector< N, T >& vec ) const
	{
		return Norm_1::norm1_impl< N, T >().operator()( vec );
	}

protected:
	/**
	 * @ingroup math
	 * Internal functor for recursive implementation.
	 *
	 * @tparam N dimension of vector
	 * @tparam T builtin-type of vector
	 */
	template< std::size_t N, typename T >
	struct norm1_impl
	{
		template< typename VecType >
		T operator() ( const VecType& vec, const T sqsum ) const
		{
			const T sq = sqsum + vec[ N-1 ] * vec[ N-1 ];
			return norm1_impl< N-1, T >().operator()( vec, sq );
		}
	};

	/**
	 * @ingroup math
	 * Partial specialization of functor for final element.
	 *
	 * @tparam N dimension of vector
	 * @tparam T builtin-type of vector
	 */
	template< typename T >
	struct norm1_impl< 0, T >
	{
		template< typename VecType >
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
 * 1st norm of a vector.
 */
struct Norm_2
{
public:
	/**
	 * @ingroup math
	 * Calculates the length of a Math::Vector.
	 *
	 * @tparam N dimension of vector
	 * @tparam T builtin-type of vector
	 * @param vec the vector 
	 * @return norm of the N-dimensional vector
	 */
	template< std::size_t N, typename T >
	T operator() ( const Math::Vector< N, T >& vec ) const
	{
		return std::sqrt( Norm_2().operator()( vec ) );
	}
};


} } } // namespace Ubitrack::Math::Functors

#endif //__H__VECTOR_FUNCTORS__

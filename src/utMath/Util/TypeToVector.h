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
 * Several casts of measurements types to vector representation
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */


#ifndef __UBITRACK_MATH_UTIL_CAST_TO_VECTOR_H_INCLUDED__
#define __UBITRACK_MATH_UTIL_CAST_TO_VECTOR_H_INCLUDED__

// Ubitrack
#include "type_traits.h"
#include "../../utUtil/Exception.h"

namespace Ubitrack { namespace Math { 

// some forward declarations for faster compiling
// class Pose;
// class Quaternion;
template< typename T, std::size_t N > class Vector;



namespace Util {

/**
 * This type_trait determines the length of
 * a vector that can represent this measurement type
 * fully and adequately.
 */
template< typename T > struct length;

template< >
struct length< float >
{
	typedef std::size_t size_type;
	static const std::size_t value = 1;
};

template< >
struct length< double >
{
	typedef std::size_t size_type;
	static const std::size_t value = 1;
};

template< >
struct length< long double >
{
	typedef std::size_t size_type;
	static const std::size_t value = 1;
};

template< typename T, std::size_t N >
struct length< Math::Vector< T, N > >
{
	typedef std::size_t size_type;
	static const std::size_t value = N;
};

template< typename T, std::size_t N >
struct length< Math::ErrorVector< T, N > >
{
	typedef std::size_t size_type;
	static const std::size_t value = N;
};

template< typename T >
struct length< Math::Scalar< T > >
{
	typedef std::size_t size_type;
	static const std::size_t value = 1;
};

template< >
struct length< Math::Quaternion >
{
	typedef std::size_t size_type;
	static const std::size_t value = 4;
};

template< >
struct length< Math::Pose >
{
	typedef std::size_t size_type;
	static const std::size_t value = 7;
};

template< >
struct length< Math::ErrorPose >
{
	typedef std::size_t size_type;
	static const std::size_t value = 7;
};

/**
 * @brief TypeToVector is a functor to cast various measurement types
 * of Ubitrack into a vector representation.
 *
 * @tparam input_type names the type of the measurement that should be casted.
 */
template< typename input_type >
struct TypeToVector
{
	/** datatype that should be casted into a vector */
	typedef input_type value_type;
	
	/** type of the value that signs the dimension of the vector. */
	typedef std::size_t size_type;
	
	/** precision of the type that will be casted, usually \c double or \c float */
	typedef typename Util::precision< value_type >::type precision_type;
	
	/** dimension of the resulting vector */
	static const size_type size = Util::length< value_type >::value;
	
	/** defines the type to which will be casted */
	typedef typename Math::Vector< precision_type, size > result_type;
	
	/** introduce an alias as the return_type is the same as the result_type */
	typedef result_type return_type;
	
	/** this binary bracket operator should be overload outside the strcut */
	void operator() ( const value_type& value, result_type& result ) const
	{
		UBITRACK_STATIC_ASSERT( ( Ubitrack::Util::is_same< value_type, result_type >::value ), YOU_NEED_TO_DEFINE_AN_OWN_CONVERSION_OPERATION );
		result = value;
		
	}

	/** this unary bracket operator is provided per default, as it just
	calls the other operator that includes the cast operations*/
	result_type operator() ( const value_type& value ) const
	{
		result_type result;
		operator()( value, result );
		return result;
	}
};

/// @internal specialization of binary bracket operator for Quaternion type
template<>
void TypeToVector< Math::Quaternion >::operator() ( const Math::Quaternion &value, result_type &rhs ) const
{
	rhs[ 0 ] = value.x();
	rhs[ 1 ] = value.y();
	rhs[ 2 ] = value.z();
	rhs[ 3 ] = value.w();
}

/// @internal specialization of unary bracket operator for Quaternion type
template<>
TypeToVector< Math::Quaternion >::result_type TypeToVector< Math::Quaternion >::operator() ( const Math::Quaternion &value ) const
{
	return result_type ( value.x(), value.y(), value.z(), value.w() );
}

/// @internal specialization of binary bracket operator for Pose type
template<>
void TypeToVector< Math::Pose >::operator() ( const Math::Pose &value, result_type &rhs ) const
{
	rhs[ 0 ] = value.translation()[ 0 ];
	rhs[ 1 ] = value.translation()[ 1 ];
	rhs[ 2 ] = value.translation()[ 2 ];
	
	rhs[ 3 ] = value.rotation().x();
	rhs[ 4 ] = value.rotation().y();
	rhs[ 5 ] = value.rotation().z();
	rhs[ 6 ] = value.rotation().w();
}

/** convenience function to support lazy use of cast operations as a binary function. The values of the measurement are set to the second parameter, a reference to the vector. */
template< typename FromType, typename ToVector >
inline void castToVector( const FromType& value, ToVector& result )
{
	UBITRACK_STATIC_ASSERT( ( Ubitrack::Util::is_same< ToVector, typename TypeToVector< FromType >::result_type >::value ), CAST_TO_DESIRED_TYPE_NOT_POSSIBLE );
	TypeToVector< FromType >()( value, result );
}

/** convenience function to support lazy use of cast operations as a unary function. The functions returns a vector containing the values of the measurement. */
template<typename ToVector, typename FromType >
inline ToVector castToVector( const FromType& value )
{
	UBITRACK_STATIC_ASSERT( ( Ubitrack::Util::is_same< ToVector, typename TypeToVector< FromType >::return_type >::value ), CAST_TO_DESIRED_TYPE_NOT_POSSIBLE );
	return TypeToVector< FromType >()( value );
}

} } } //namespace Ubitrack::Math::Util

#endif // __UBITRACK_MATH_UTIL_CAST_TO_VECTOR_H_INCLUDED__
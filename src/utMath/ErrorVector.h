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
 * @file
 * A vector with an associated covariance matrix.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_MATH_ERRORVECTOR_H_INCLUDED__
#define __UBITRACK_MATH_ERRORVECTOR_H_INCLUDED__


#include "Vector.h"
#include "Matrix.h"

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/blas/blas3.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

namespace Ubitrack { namespace Math {

/**
 * This class stores an N-vector with an associated covariance matrix 
 * and provides methods for transformations, etc.
 *
 * @tparam T type of vector/matrix elements
 * @tparam N size of the vector
 */
template< typename T, std::size_t N >
struct ErrorVector
{
	typedef T value_type;
	
	/** vector contents */
	Math::Vector< T, N > value;

	/** covariance matrix of \c value */
	Math::Matrix< T, N, N > covariance;
	
	/** default constructor */
	ErrorVector()
	{}

	/** construct from value and covariance */
	template< class VT, class MT >
	ErrorVector( const VT& _value, const MT& _covariance )
		: value( _value )
		, covariance( _covariance )
	{}


	/** compute RMS value (in this case: square-root of trace of covariance matrix) */
	T getRMS( void ) 
	{
		T trace = 0;
		for ( std::size_t i = 0; i < N; i++ ) 
		{
			trace += covariance( i,i );
		}
		
		return sqrt ( trace );
	}
	
	
	/**
	 * (un-)serialization helper function
	 */
	template< class Archive > 
	void serialize( Archive& ar, const unsigned int version )
	{
		ar & value;
		ar & covariance;
	}
};


/// @internal stream output operator
template< typename T, std::size_t N >
std::ostream& operator<<( std::ostream& s, const ErrorVector< T, N >& v )
{
	s << v.value << std::endl << v.covariance;
	return s;
}


} } // namespace Ubitrack::Math

#endif //__UBITRACK_MATH_ERRORVECTOR_H_INCLUDED__

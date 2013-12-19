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
 * @ingroup Math
 * @file
 * vector class to store intermediate results in.
 * specialized for fixed-size vectors with optimized stack-storage.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_DETAIL_RESULTVECTOR_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_DETAIL_RESULTVECTOR_H_INCLUDED__

#include <utMath/Vector.h>

namespace Ubitrack { namespace Math { namespace Optimization { namespace Function { namespace Detail {

/**
 * statically sized, stack-allocated vector for intermediate results
 */
template< std::size_t Size >
class ResultVector
	: public Math::Vector< double, Size >
{
public:
	ResultVector( std::size_t s )
	{
		assert( s == Size );
	}

	template< class AE > 
	ResultVector( const boost::numeric::ublas::vector_expression< AE >& e )
		: Math::Vector< double, Size >( e )
	{}
	
	template< class AE > 
	Math::Vector< double, Size >& operator=( const boost::numeric::ublas::vector_expression< AE >& e )
	{
		Math::Vector< double, Size >::operator=( e );
		return *this;
	}
};

/** dynamically sized, heap-allocated vector for intermediate results */
template<>
class ResultVector< 0U >
	: public Math::Vector< double >
{
public:
	ResultVector( std::size_t s )
		: Math::Vector< double >( s )
	{}

	template< class AE > 
	ResultVector( const boost::numeric::ublas::vector_expression< AE >& e )
		: Math::Vector< double >( e )
	{}
	
	template< class AE > 
	Math::Vector< double >& operator=( const boost::numeric::ublas::vector_expression< AE >& e )
	{
		Math::Vector< double >::operator=( e );
		return *this;
	}
};

}}}}} // namespace Ubitrack::Math::Optimization::Function::Detail

#endif

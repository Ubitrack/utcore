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
 * matrix class to store intermediate results in.
 * specialized for fixed-size matrices with optimized stack-storage.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_DETAIL_RESULTMATRIX_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_DETAIL_RESULTMATRIX_H_INCLUDED__

#include <utMath/Matrix.h>

namespace Ubitrack { namespace Math { namespace Optimization { namespace Function { namespace Detail {

/**
 * statically sized, stack-allocated matrix for intermediate results
 */
template< std::size_t TotalSize, std::size_t Size1, std::size_t Size2 >
class ResultMatrix
	: public Math::Matrix< double, Size1, Size2 >
{
public:
	ResultMatrix( std::size_t s1, std::size_t s2 )
	{
		assert( Size1 == s1 && Size2 == s2 );
	}
	
	template< class AE > 
	Math::Matrix< double, Size1, Size2 >& operator=( const boost::numeric::ublas::matrix_expression< AE >& e )
	{
		Math::Matrix< double, Size1, Size2 >::operator=( e );
		return *this;
	}
};

/** dynamically sized, heap-allocated matrix for intermediate results */
template< std::size_t Size1, std::size_t Size2 >
class ResultMatrix< Size2, 0U, Size1 >
	: public Math::Matrix< double >
{
public:
	ResultMatrix( std::size_t s1, std::size_t s2 )
		: Math::Matrix< double >( s1, s2 )
	{}
	
	template< class AE > 
	Math::Matrix< double >& operator=( const boost::numeric::ublas::matrix_expression< AE >& e )
	{
		Math::Matrix< double >::operator=( e );
		return *this;
	}
};

}}}}} // namespace Ubitrack::Math::Optimization::Function::Detail

#endif

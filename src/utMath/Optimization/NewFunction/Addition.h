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
 * adds two values
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_ADDITION_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_ADDITION_H_INCLUDED__
 
#include <utMath/Optimization/NewFunction/MultiVariateFunction.h>
 
namespace Ubitrack { namespace Math { namespace Optimization { namespace Function {

/**
 * Function that adds two vectors.
 */
template< std::size_t N >
class Addition
	: public MultiVariateFunction< Addition< N >, N >
{
public:

	template< class DestinationVector, class Param1, class Param2 >
	void evaluate( DestinationVector& result, const Param1& p1, const Param2& p2 ) const
	{
		result = p1 + p2;
	}
		
	template< class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void multiplyJacobian1( const LeftHand& l, DestinationMatrix& j, const Param1&, const Param2& ) const
	{
		j = l;
	}
	
	template< class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void multiplyJacobian2( const LeftHand& l, DestinationMatrix& j, const Param1&, const Param2& ) const
	{
		j = l;
	}
};

}}}} // namespace Ubitrack::Math::Optimization::Function

#endif

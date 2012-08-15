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
 * stores an intermediate result
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_STOREINTERMEDIATE_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_STOREINTERMEDIATE_H_INCLUDED__
 
#include <utMath/NewFunction/MultiVariateFunction.h>
 
namespace Ubitrack { namespace Math { namespace Function {

/**
 * Stores an intermediate result in a vector
 */
template< unsigned M, class T = double >
class StoreIntermediate
	: public MultiVariateFunction< StoreIntermediate< M, T >, M >
{
public:
	/** 
	 * Construct from matrix. 
	 * Note: matrix reference must be valid throughout the lifetime of the object! 
	 */
	StoreIntermediate( Math::Vector< M, T >& _vector )
		: rVector( _vector )
	{}

	template< class DestinationVector, class Param1 >
	void evaluate( DestinationVector& result, const Param1& p1 ) const
	{
		result = p1;
		rVector = p1;
	}
		
	template< class LeftHand, class DestinationMatrix, class Param1 >
	void multiplyJacobian1( const LeftHand& l, DestinationMatrix& j, const Param1& ) const
	{
		j = l;
	}
	
	mutable Math::Vector< M, T >& rVector;
};

} } } // namespace Ubitrack::Math::Function

#endif

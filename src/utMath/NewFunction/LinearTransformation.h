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
 * linear transformation of a vector, assuming the transformation matrix is constant
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_LINEARTRANSFORMATION_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_LINEARTRANSFORMATION_H_INCLUDED__
 
#include <utMath/NewFunction/MultiVariateFunction.h>
 
namespace Ubitrack { namespace Math { namespace Function {

/**
 * Function that multiplies a vector with a constant matrix
 */
template< unsigned M, unsigned N, class T = double >
class LinearTransformation
	: public MultiVariateFunction< LinearTransformation< M, N, T >, M >
{
public:
	/** 
	 * Construct from matrix. 
	 * Note: matrix reference must be valid throughout the lifetime of the object! 
	 */
	LinearTransformation( const Math::Matrix< M, N, T >& _matrix )
		: rMatrix( _matrix )
	{}

	template< class DestinationVector, class Param1 >
	void evaluate( DestinationVector& result, const Param1& p1 ) const
	{
		result = boost::numeric::ublas::prod( rMatrix, p1 );
	}
		
	template< class LeftHand, class DestinationMatrix, class Param1 >
	void multiplyJacobian1( const LeftHand& l, DestinationMatrix& j, const Param1& ) const
	{
		j = boost::numeric::ublas::prod( l, rMatrix );
	}
	
	const Math::Matrix< M, N, T >& rMatrix;
};

} } } // namespace Ubitrack::Math::Function

#endif

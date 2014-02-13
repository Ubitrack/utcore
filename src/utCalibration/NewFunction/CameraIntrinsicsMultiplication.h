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
 * @ingroup calibration
 * @file
 * Multiplies a 2-vector with a camera intrinsics matrix given as a 5-vector
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_CAMERAINTRINSICSMULTIPLICATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_CAMERAINTRINSICSMULTIPLICATION_H_INCLUDED__

#include <utMath/Optimization/NewFunction/MultiVariateFunction.h>
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Applies an intrinsic camera matrix (given as a 5-vector) to an already dehomogenized 2-vector.
 * We assume that the lower right matrix element is -1. See src/Vision/Readme.txt why.
 */
class CameraIntrinsicsMultiplication
	: public Math::Optimization::Function::MultiVariateFunction< CameraIntrinsicsMultiplication, 2 >
{
public:

	template< class DestinationVector, class Param1, class Param2 >
	void evaluate( DestinationVector& result, const Param1& intr, const Param2& point ) const
	{
		result( 0 ) = - ( intr( 0 ) * point( 0 ) + intr( 1 ) * point( 1 ) + intr( 2 ) );
		result( 1 ) = - (                          intr( 3 ) * point( 1 ) + intr( 4 ) );
	}

	/** calculate jacobian wrt the matrix */
	template< class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void multiplyJacobian1( const LeftHand& l, DestinationMatrix& j, const Param1& intr, const Param2& point ) const
	{
		const std::size_t n_rows( l.size1() );
		for ( std::size_t r ( 0 ); r < n_rows; r++ )
		{
			j( r, 0 ) = -point( 0 ) * l( r, 0 );
			j( r, 1 ) = -point( 1 ) * l( r, 0 );
			j( r, 2 ) = -             l( r, 0 );
			j( r, 3 ) = -point( 1 ) * l( r, 1 );
			j( r, 4 ) = -             l( r, 1 );
		}

		/*J( 0, 0 ) = -point( 0 );
		J( 0, 1 ) = -point( 1 );
		J( 0, 2 ) = -1;
		J( 0, 3 ) = 0;
		J( 0, 4 ) = 0;
		J( 1, 0 ) = 0;
		J( 1, 1 ) = 0;
		J( 1, 2 ) = 0;
		J( 1, 3 ) = -point( 1 );
		J( 1, 4 ) = -1;*/
	}

	/** calculate jacobian wrt the point */
	template< class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void multiplyJacobian2( const LeftHand& l, DestinationMatrix& j, const Param1& intr, const Param2& point ) const
	{
		const std::size_t n_rows( l.size1() );
		for ( std::size_t r ( 0 ); r < n_rows; r++ )
		{
			j( r, 0 ) = -  intr( 0 ) * l( r, 0 );
			j( r, 1 ) = -( intr( 1 ) * l( r, 0 ) + intr( 3 ) * l( r, 1 ) );
		}

		/*K22( 0, 0 ) = -input( iP ); K22( 0, 1 ) = -input( iP + 1 );
		K22( 1, 0 ) = 0;            K22( 1, 1 ) = -input( iP + 3 );*/
	}
};

}}} // namespace Ubitrack::Calibration::Function

#endif

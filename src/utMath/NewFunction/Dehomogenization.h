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
 * functions for division by the last element
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_DEHOMOGENIZATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_DEHOMOGENIZATION_H_INCLUDED__

#include <utMath/NewFunction/MultiVariateFunction.h>
 
namespace Ubitrack { namespace Math { namespace Function {

/**
 * Function that dehomogenizes a vector by dividing through the last element and then dropping it.
 * Given an N-vector it returns an (N-1)-Vector
 */
template< std::size_t N >
struct Dehomogenization
	: public MultiVariateFunction< Dehomogenization< N >, N-1 >
{
	/*
	 * @param result an N-1-vector
	 * @param input an N-vector
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		typename VT1::value_type f( 1 / input( N - 1 ) );
		for ( std::size_t i( 0 ); i < N - 1; i++ )
			result( i ) = input( i ) * f;
	}
	
	template< class LeftHand, class DestinationMatrix, class Param1 >
	void multiplyJacobian1( const LeftHand& l, DestinationMatrix& j, const Param1& input ) const
	{
		typename LeftHand::value_type tz = 1 / input( N-1 );
		typename LeftHand::value_type tz2 = tz * tz;
		for ( std::size_t r = 0; r < l.size1(); r++ )
		{
			for ( std::size_t c = 0; c < (N-1); c++ )
				j( r, c ) = tz * l( r, c );
			j( r, N-1 ) = -boost::numeric::ublas::inner_prod( boost::numeric::ublas::row( l, r ), boost::numeric::ublas::subrange( input, 0, N-1 ) ) * tz2;
			
		}

		// TODO: check for division by zero
		/*typename MT::value_type tz = 1 / input( N-1 );
		boost::numeric::ublas::subrange( J, 0, N-1, 0, N-1 ) = 
			boost::numeric::ublas::identity_matrix< typename MT::value_type >( N-1 ) * tz;
		tz = tz * tz;
		for ( unsigned i = 0; i < N-1; i++ )
			J( i, N-1 ) = -input( i ) * tz;*/
	}
};

} } } // namespace Ubitrack::Calibration::Function

#endif

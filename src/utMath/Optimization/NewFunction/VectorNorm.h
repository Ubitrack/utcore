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
 * Defines a function that computes the norm of a vector
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#include <utMath/Vector.h>

namespace Ubitrack { namespace Math { namespace Optimization { namespace Function {

/**
 * A function that computes the norm of a vector
 */
template< unsigned M >
struct VectorNorm
	: public MultiVariateFunction< VectorNorm< M >, 1 >
{
	template< class DestinationVector, class Param1 >
	void evaluate( DestinationVector& result, const Param1& p1 ) const
	{
		result( 0 ) = boost::numeric::ublas::norm_2( p1 );
	}
		
	template< class LeftHand, class DestinationMatrix, class Param1 >
	void multiplyJacobian1( const LeftHand& l, DestinationMatrix& j, const Param1& p1 ) const
	{
		typename Param1::value_type f = l( 0, 0 ) / boost::numeric::ublas::norm_2( p1 );
		for ( unsigned i = 0; i < M; i++ )
			j( 0, i ) = f * p1( i );
	}
};

} } } // namespace Ubitrack::Calibration::Normalize

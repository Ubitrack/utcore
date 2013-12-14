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
 * Function that inverts a rotation velocity dR, given dR and R.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 


#ifndef __UBITRACK_TRACKING_FUNCTION_INVERTROTATIONVELOCITY_H_INCLUDED__
#define __UBITRACK_TRACKING_FUNCTION_INVERTROTATIONVELOCITY_H_INCLUDED__

#include <utMath/Vector.h>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace Ubitrack { namespace Tracking { namespace Function {

/**
 * Function that inverts a rotation velocity dR, given dR and R.
 * The resulting dR' is computed as dR' = R dR R*
 */
class InvertRotationVelocity
{
public:
	unsigned size() const
	{ return 3; }
	
	/**
	 * Updates the vector and computes the jacobian.
	 * The input vector must be of size 3 and contain (q, dq), where q is the absolute
	 * rotation value quaternion and dq the rotation velocity as a 3-vector
	 * @param result vector to put the result in
	 * @param input input vector 
	 * @param jacobian matrix to put the jacobian into
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		typedef typename VT1::value_type T;
		namespace ublas = boost::numeric::ublas;

		// invert rotation velocity
		ublas::vector< T > invertedVelocity( -ublas::subrange( input, 4, 7 ) );

		// compute results and jacobians
		ublas::matrix_range< MT > jRot( jacobian, ublas::range( 0, 3 ), ublas::range( 0, 4 ) );
		ublas::matrix_range< MT > jVel( jacobian, ublas::range( 0, 3 ), ublas::range( 4, 7 ) );
		Math::Function::QuaternionVectorRotation().evaluateWithJacobian( result, ublas::subrange( input, 0, 4 ), 
			invertedVelocity, jRot, jVel );

		// adjust jacobian for inversion of q and dq
		jVel *= T( -1 );
	}

};

} } } // namespace Ubitrack::Tracking::Function

#endif

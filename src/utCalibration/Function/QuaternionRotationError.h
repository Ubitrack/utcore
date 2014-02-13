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
 * function that rotates a vector around a quaternion+error
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_QUATERNIONROTATIONERROR_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_QUATERNIONROTATIONERROR_H_INCLUDED__
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Function that rotates a vector \c v around a quaternion \c q with associated error qe, 
 * producing v2 = (q * qe * v * qe' * q') and/or computes the jacobian of the product 
 * wrt. qe = (x, y, z, 1), assuming that qe = (0,0,0,1).
 */
template< class VType >
class QuaternionRotationError
{
public:
	/** 
	 * constructor.
	 * @param v reference to the vector to rotate (must stay constant during the lifetime of the function object
	 */
	QuaternionRotationError( const Math::Vector< VType, 3 >& v )
		: m_v( v )
	{}
	
	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 3; }
	
	/**
	 * This jacobian is only needed for error propagation, so only the jacobian() method is necessary.
	 * 
	 * @param input a 4-vector containing the quaternion q = (x, y, z, w)
	 * @param J a 3x3 matrix where the resulting jacobian (wrt. error) is stored
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& j ) const
	{
		VType t203 = input( 1 );
		VType t201 = m_v( 1 );
		VType t200 = input( 2 );
		VType t204 = m_v( 2 );
		VType t208 = input( 3 );
		VType t199 = input( 0 );
		VType t215 = m_v( 0 );
		VType t218 = t203 * t203;
		VType t219 = 2 * t218;
		VType t220 = t200 * t200;
		VType t221 = 2 * t220;
		VType t222 =  -1 + t219 + t221;
		VType t210 = t200 * t204;
		VType t232 = t199 * t199;
		VType t233 = 2 * t232;
		VType t234 =  -1 + t221 + t233;
		VType t241 = t199 * t215;
		VType t209 = t201 * t203;
		j( 0, 0 ) = 4 * ( t199 * ( t200 * t201 - t203 * t204 ) + t208 * ( t209 + t210 ) );
		j( 0, 1 ) = -2 * ( 2 * t199 * t200 * t215 + 2 * t203 * t208 * t215 + t204 * t222 );
		j( 0, 2 ) = 4 * t199 * t203 * t215 - 4 * t200 * t208 * t215 + 2 * t201 * t222;
		j( 1, 0 ) = 4 * t200 * t201 * t203 - 4 * t199 * t201 * t208 + 2 * t204 * t234;
		j( 1, 1 ) = 4 * ( t203 * ( t199 * t204 - t200 * t215 ) + t208 * ( t210 + t241 ) );
		j( 1, 2 ) = -2 * ( 2 * t201 * ( t199 * t203 + t200 * t208 ) + t215 * t234 );
		j( 2, 0 ) = -2 * ( 2 * t204 * ( t200 * t203 + t199 * t208 ) + t201 * ( -1 + t219 + t233 ) );
		j( 2, 1 ) = 4 * t204 * ( t199 * t200 - t203 * t208 ) + t215 * ( -2 + 4 * t218 + 4 * t232 );
		j( 2, 2 ) = 4 * ( t200 * (  - ( t199 * t201 ) + t203 * t215 ) + t208 * ( t209 + t241 ) );	
	}
	
protected:
	const Math::Vector< VType, 3 >& m_v;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

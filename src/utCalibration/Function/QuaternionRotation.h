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
 * functions rotations around a quaternion
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_QUATERNIONROTATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_QUATERNIONROTATION_H_INCLUDED__
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Function that rotates a vector \c v around a quaternion \c q, producing v2 = (q * v * q') 
 * and/or computes the jacobian of the product wrt. q = (x, y, z, w).
 */
template< class VType >
class QuaternionRotation
{
public:
	/** 
	 * constructor.
	 * @param v reference to the vector to rotate (must stay constant during the lifetime of the function object
	 */
	QuaternionRotation( const Math::Vector< 3, VType >& v )
		: m_v( v )
	{}
	
	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 3; }

	/**
	 * @param result a 3-vector containing the rotated vector
	 * @param input a 4-vector containing the quaternion q = (x, y, z, w)
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		result = Math::Quaternion::fromVector( input ) * m_v;
	}
		
	/**
	 * @param result a 3-vector containing the rotated vector
	 * @param input a 4-vector containing the quaternion q = (x, y, z, w)
	 * @param J a 3x4 matrix where the resulting jacobian is stored
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input a 4-vector containing the quaternion q = (x, y, z, w)
	 * @param J a 3x4 matrix where the resulting jacobian is stored
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		VType t4  = VType( input( 0 ) * m_v( 0 ) + input( 1 ) * m_v( 1 ) + input( 2 ) * m_v( 2 ) );
		VType t8  = VType( input( 3 ) * m_v( 2 ) + input( 0 ) * m_v( 1 ) - input( 1 ) * m_v( 0 ) );
		VType t12 = VType( input( 0 ) * m_v( 2 ) - input( 3 ) * m_v( 1 ) - input( 2 ) * m_v( 0 ) );
		VType t16 = VType( input( 3 ) * m_v( 0 ) + input( 1 ) * m_v( 2 ) - input( 2 ) * m_v( 1 ) );
		J( 0, 0 ) = 2 * t4;
		J( 0, 1 ) = 2 * t8;
		J( 0, 2 ) = 2 * t12;
		J( 0, 3 ) = 2 * t16;
		J( 1, 0 ) = -2 * t8;
		J( 1, 1 ) = 2 * t4;
		J( 1, 2 ) = 2 * t16;
		J( 1, 3 ) = -2 * t12;
		J( 2, 0 ) = -2 * t12;
		J( 2, 1 ) = -2 * t16;
		J( 2, 2 ) = 2 * t4;
		J( 2, 3 ) = 2 * t8;
	}
	
protected:
	const Math::Vector< 3, VType >& m_v;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

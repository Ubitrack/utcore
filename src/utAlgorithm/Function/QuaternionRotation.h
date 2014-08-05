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

#ifndef __UBITRACK_ALGORITHM_FUNCTION_QUATERNIONROTATION_H_INCLUDED__
#define __UBITRACK_ALGORITHM_FUNCTION_QUATERNIONROTATION_H_INCLUDED__
 
namespace Ubitrack { namespace Algorithm { namespace Function {

/**
 * Function that rotates a vector \c v around a quaternion \c q, producing v2 = (q * v * q') 
 * and/or computes the jacobian of the product wrt. q = (x, y, z, w).
 */
template< class VType >
class QuaternionRotation
{
protected:
	const Math::Vector< VType, 3 >& m_v;
	
public:
	/** 
	 * constructor.
	 * @param v reference to the vector to rotate (must stay constant during the lifetime of the function object
	 */
	QuaternionRotation( const Math::Vector< VType, 3 >& v )
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
	
		/*
		%% short matlab symbolic expression example code to 
		%% illustrate derivation of the jacobian for the vector rotation
		syms qx qy qz qw x y z
		mat3x3 = [ (qw*qw + qx*qx - qy*qy - qz*qz), (2*qx*qy - 2*qw*qz) , (2*qx*qz + 2*qw*qy);
				(2*qx*qy + 2*qw*qz), (qw*qw - qx*qx + qy*qy - qz*qz ), (2*qy*qz - 2*qw*qx);
				(2*qx*qz - 2*qw*qy) , (2*qy*qz + 2*qw*qx) , qw*qw - qx*qx - qy*qy + qz*qz]
		f = mat3x3 * [x; y; z]
		jacobian( (f) , [qx, qy, qz] ) %qw is skipped
		*/
		const VType qx = input[ 0 ];
		const VType qy = input[ 1 ];
		const VType qz = input[ 2 ];
		const VType qw = input[ 3 ];
		const VType x = m_v[ 0 ];
		const VType y = m_v[ 1 ];
		const VType z = m_v[ 2 ];
		
		const VType t2 = qx*y*2;
		const VType t3 = qw*z*2;
		const VType t4 = qx*x*2;
		const VType t5 = qy*y*2;
		const VType t6 = qz*z*2;
		const VType t7 = t4+t5+t6;
		const VType t8 = qw*x*2;
		const VType t9 = qy*z*2;
		const VType t15 = qz*y*2;
		const VType t10 = t8+t9-t15;
		const VType t11 = qx*z*2;
		const VType t12 = qz*x*2;
		const VType t13 = qw*y*2;
		const VType t14 = -t11+t12+t13;
		const VType t16 = qy*x*2;
		J( 0, 0 ) = t7;
		J( 0, 1 ) = t2+t3-qy*x*2.0;
		J( 0, 2 ) = t11-qz*x*2.0-qw*y*2.0;
		J( 0, 3 ) = t10;
		J( 1, 0 ) = -t2-t3+t16;
		J( 1, 1 ) = t7;
		J( 1, 2 ) = t10;
		J( 1, 3 ) = t14;
		J( 2, 0 ) = t14;
		J( 2, 1 ) = -t8-t9+t15;
		J( 2, 2 ) = t7;
		J( 2, 3 ) = t2+t3-t16;
		
		// old code, left for comparison
		// VType t4  = VType( input( 0 ) * m_v( 0 ) + input( 1 ) * m_v( 1 ) + input( 2 ) * m_v( 2 ) );
		// VType t8  = VType( input( 3 ) * m_v( 2 ) + input( 0 ) * m_v( 1 ) - input( 1 ) * m_v( 0 ) );
		// VType t12 = VType( input( 0 ) * m_v( 2 ) - input( 3 ) * m_v( 1 ) - input( 2 ) * m_v( 0 ) );
		// VType t16 = VType( input( 3 ) * m_v( 0 ) + input( 1 ) * m_v( 2 ) - input( 2 ) * m_v( 1 ) );
		// J( 0, 0 ) = 2 * t4;
		// J( 0, 1 ) = 2 * t8;
		// J( 0, 2 ) = 2 * t12;
		// J( 0, 3 ) = 2 * t16;
		// J( 1, 0 ) = -2 * t8;
		// J( 1, 1 ) = 2 * t4;
		// J( 1, 2 ) = 2 * t16;
		// J( 1, 3 ) = -2 * t12;
		// J( 2, 0 ) = -2 * t12;
		// J( 2, 1 ) = -2 * t16;
		// J( 2, 2 ) = 2 * t4;
		// J( 2, 3 ) = 2 * t8;
	}
};

} } } // namespace Ubitrack::Algorithm::Function

#endif //__UBITRACK_ALGORITHM_FUNCTION_QUATERNIONROTATION_H_INCLUDED__

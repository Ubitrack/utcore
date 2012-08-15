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
 * rotates a vector around a 3-element lie group rotation
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_LIEROTATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_LIEROTATION_H_INCLUDED__
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Function that rotates a vector \c v around a quaternion \c q, producing v2 = (q * v * q'). The quaternion q is 
 * represented by its 3-element logarithm (r_x, r_y, r_z).
 * The function class computes the jacobian of the product wrt. r = (r_x, r_y, r_z).
 */
template< class VType >
class LieRotation
{
public:
	/** 
	 * constructor.
	 * @param v reference to the vector to rotate (must stay constant during the lifetime of the function object
	 */
	LieRotation( const Math::Vector< 3, VType >& v )
		: m_v( v )
	{}
	
	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 3; }

	/**
	 * @param result a 3-vector containing the rotated vector
	 * @param input a 3-vector containing the rotation r = (r_x, r_y, r_z)
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		result = Math::Quaternion::fromLogarithm( input ) * m_v;
	}
		
	/**
	 * @param result a 3-vector containing the rotated vector
	 * @param input a 3-vector containing the rotation r = (r_x, r_y, r_z)
	 * @param J a 3x3 matrix where the resulting jacobian is stored
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input a 3-vector containing the rotation r = (r_x, r_y, r_z)
	 * @param J a 3x3 matrix where the resulting jacobian is stored
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		VType t51944 = input( 0 );
		VType t51945 = t51944 * t51944;
		VType t51946 = input( 1 );
		VType t51947 = t51946 * t51946;
		VType t51948 = input( 2 );
		VType t51949 = t51948 * t51948;
		VType t51950 = t51945 + t51947 + t51949;
		VType t51952 = sqrt( t51950 );
		
		if ( t51952 < VType( 1e-6 ) )
		{
			// avoid division by zero and compute limit of J with r->0
			J( 0, 0 ) = 0;
			J( 0, 1 ) = m_v( 2 );
			J( 0, 2 ) = -m_v( 1 );
			J( 1, 0 ) = -m_v( 2 );
			J( 1, 1 ) = 0;
			J( 1, 2 ) = m_v( 0 );
			J( 2, 0 ) = m_v( 1 );
			J( 2, 1 ) = -m_v( 0 );
			J( 2, 2 ) = 0;
		}
		else
		{
			VType t51953 = cos( t51952 );
			VType t51955 = m_v( 1 );
			VType t51957 = m_v( 2 );
			VType t51963 = t51946 * t51955;
			VType t51964 = t51948 * t51957;
			VType t51965 = t51963 + t51964;
			VType t51962 = t51947 + t51949;
			VType t51967 = sin( t51952 );
			VType t51971 = m_v( 0 );
			VType t51956 =  - ( t51948 * t51955 );
			VType t51958 = t51946 * t51957;
			VType t51951 = 1 / ( t51950 * t51950 );
			VType t51954 = t51944 * t51944 * t51944;
			VType t51959 = t51956 + t51958;
			VType t52004 =  - t51949;
			VType t52001 = 2 * t51971;
			VType t52002 = t51956 + t51958 + t52001;
			VType t51961 =  -1 + t51953;
			VType t52058 =  - ( t51946 * t51971 );
			VType t52015 = t51948 * t51948 * t51948;
			VType t51984 = t51946 * t51946 * t51946;
			VType t52022 = t51944 * t51971;
			VType t52023 = t51963 + t52022;
			VType t52024 =  - 2 * t51948 * t52023;
			VType t52025 = t51945 + t51947 + t52004;
			VType t52026 = t51957 * t52025;
			VType t52027 = t52024 + t52026;
			VType t52045 = ( t51944 * t51944 ) * ( t51944 * t51944 );
			VType t52042 =  - ( t51952 * t51955 * t51967 );
			VType t52083 =  - ( t51948 * t51971 );
			VType t52090 = t51948 * t51971;
			VType t52048 = t51952 * t51957 * t51967;
			VType t52059 = t51946 * t51952 * t51967 * t51971;
			VType t51988 = t51946 * t51948 * t51952 * t51955 * t51967;
			VType t52110 = t51944 * t51946 * t51948 * t51952 * t51967 * t51971;
			VType t51974 =  - ( t51949 * t51952 * t51967 * t51971 );
			VType t51992 =  -2 * t51944 * t51946 * t51971;
			VType t51993 = t51945 * t51955;
			VType t51994 =  - ( t51947 * t51955 );
			VType t51995 = t51949 * t51955;
			VType t51996 =  -2 * t51946 * t51948 * t51957;
			VType t51997 = t51992 + t51993 + t51994 + t51995 + t51996;
			VType t52139 = t51955 * t51962;
			VType t52067 = t51954 * t51971;
			VType t52191 = t51956 + t51971;
			VType t52167 = t51948 * t51955;
			VType t52198 = t51955 + t52090;
			VType t52007 = 2 * t51946 * t51948 * t51957;
			J( 0, 0 ) = t51951 * ( t51953 * t51954 * t51959 - t51961 * t51962 * t51965 + t51945 * t51965 * (  - 1 + t51953 + t51952 * t51967 ) + 
				t51944 * ( t51948 * t51952 * t51955 * t51967 - t51946 * t51952 * t51957 * t51967 + 
				t51953 * t51962 * ( t51956 + t51958 - 2 * t51971 ) + 2 * t51962 * t51971 - t51947 * t51952 * t51967 * t51971 + t51974 ) );
			J( 0, 1 ) = t51951 * ( t51944 * t51947 * t51952 * t51955 * t51967 + t51945 * t51952 * t51957 * t51967 + 
				t51944 * t51946 * t51948 * t51952 * t51957 * t51967 + t51949 * t51952 * t51957 * t51967 - 
				t51946 * t51949 * t51952 * t51967 * t51971 - t51952 * t51967 * t51971 * t51984 + t51988 + t51944 * t51997 + 
				t51953 * (  - ( t51954 * t51955 ) + t51946 * t51959 * t51962 + t51945 * t51946 * t52002 + 
				t51944 * ( t51955 * ( t51947 + t52004 ) + t52007 ) ) );
			J( 0, 2 ) = t51951 * (  - ( t51945 * t51952 * t51955 * t51967 ) - t51947 * t51952 * t51955 * t51967 + 
				t51944 * t51946 * t51948 * t51952 * t51955 * t51967 - t51946 * t51948 * t51952 * t51957 * t51967 + 
				t51944 * t51949 * t51952 * t51957 * t51967 - t51947 * t51948 * t51952 * t51967 * t51971 + 
				t51953 * (  - ( t51954 * t51957 ) + t51944 * ( 2 * t51946 * t51948 * t51955 - t51947 * t51957 + t51949 * t51957 ) + 
				t51948 * t51959 * t51962 + t51945 * t51948 * t52002 ) - t51952 * t51967 * t51971 * t52015 + t51944 * t52027 );
			J( 1, 0 ) = t51951 * ( t51944 * ( 2 * t51946 * t51961 * t51965 - t51949 * t51952 * t51955 * t51967 + 
				t51946 * t51948 * t51952 * t51957 * t51967 + t51948 * t51953 * t51962 * t51971 - t51948 * t51952 * t51967 * t51971 ) + 
				t51954 * ( t51948 * t51953 * t51971 + t52042 ) - t51953 * t51957 * t52045 - t51962 * ( t51946 * t51961 * t51971 + t52048 ) + 
				t51945 * ( t52058 - t51953 * ( t51957 * t51962 + t52058 ) + t52059 ) );
			J( 1, 1 ) = t51951 * ( 2 * t51945 * t51946 * t51955 + 2 * t51946 * t51949 * t51955 + t51945 * t51948 * t51957 - 
				t51947 * t51948 * t51957 - t51945 * t51946 * t51952 * t51955 * t51967 - t51946 * t51949 * t51952 * t51955 * t51967 + 
				t51944 * t51946 * t51952 * t51957 * t51967 + t51947 * t51948 * t51952 * t51957 * t51967 - t51944 * t51947 * t51971 + 
				t51944 * t51949 * t51971 + t51944 * t51947 * t51952 * t51967 * t51971 - t51946 * t51948 * t51952 * t51967 * t51971 + 
				t51957 * t52015 + t52067 - t51953 * ( t51954 * ( t51958 + t51971 ) + 
				t51944 * ( t51946 * t51957 * t51962 + (  - t51947 + t51949 ) * t51971 ) + 
				t51945 * ( t51964 + t51946 * ( 2 * t51955 + t52083 ) ) - 
				t51948 * ( t51947 * t51957 - t51949 * t51957 + t51971 * t51984 + t51946 * t51948 * (  - 2 * t51955 + t52090 ) ) ) );
			J( 1, 2 ) = t51951 * (  - ( t51945 * t51948 * t51952 * t51955 * t51967 ) + t51944 * t51948 * t51952 * t51957 * t51967 + 
				t51946 * t51949 * t51952 * t51957 * t51967 + t51945 * t51952 * t51967 * t51971 + t51947 * t51952 * t51967 * t51971 - 
				t51952 * t51955 * t51967 * t52015 + 
				t51953 * (  - ( t51957 * ( t51945 * t51946 - t51946 * t51949 + t51948 * t51954 + t51944 * t51948 * t51962 + t51984 ) ) + 
				t51948 * ( 2 * t51947 * t51955 + t51971 * ( 2 * t51944 * t51946 + ( t51945 + t51947 ) * t51948 + t52015 ) ) ) + 
				t51946 * t52027 + t52110 );
			J( 2, 0 ) = t51951 * (  - ( t51962 * ( t51948 * t51961 * t51971 + t52042 ) ) + t51953 * t51955 * t52045 - 
				t51954 * ( t51946 * t51953 * t51971 + t52048 ) + 
				t51944 * ( 2 * t51948 * t51961 * t51965 - t51947 * t51952 * t51957 * t51967 - t51946 * t51953 * t51962 * t51971 + t51988 + t52059 ) + 
				t51945 * ( t51948 * t51952 * t51967 * t51971 + t52083 + t51953 * ( t52090 + t52139 ) ) );
			J( 2, 1 ) = t51951 * (  - ( t51944 * t51946 * t51952 * t51955 * t51967 ) + t51947 * t51948 * t51952 * t51955 * t51967 - 
				t51945 * t51946 * t51952 * t51957 * t51967 - t51945 * t51952 * t51967 * t51971 + t51974 - t51952 * t51957 * t51967 * t51984 + 
				t51948 * t51997 + t52110 + t51953 * ( t51946 * t51954 * t51955 + 2 * t51946 * t51949 * t51957 - t51955 * t52015 + 
				t51947 * t51948 * ( t51955 + t52083 ) + t51944 * t51946 * ( 2 * t51948 * t51971 + t52139 ) - 
				t51945 * ( t51947 * t51971 + t52167 ) - t51971 * ( t51946 * t51946 ) * ( t51946 * t51946 ) ) );
			J( 2, 2 ) = t51951 * ( t51945 * t51946 * t51955 - t51946 * t51949 * t51955 + 2 * t51945 * t51948 * t51957 + 
				2 * t51947 * t51948 * t51957 - t51944 * t51948 * t51952 * t51955 * t51967 + t51946 * t51949 * t51952 * t51955 * t51967 - 
				t51945 * t51948 * t51952 * t51957 * t51967 - t51947 * t51948 * t51952 * t51957 * t51967 + t51944 * t51947 * t51971 - 
				t51944 * t51949 * t51971 + t51946 * t51948 * t51952 * t51967 * t51971 + t51944 * t51949 * t51952 * t51967 * t51971 + 
				t51955 * t51984 + t52067 - t51953 * ( t51954 * t52191 + t51944 * (  - ( t51949 * ( t51971 + t52167 ) ) + t51947 * t52191 ) + 
				t51945 * ( 2 * t51948 * t51957 + t51946 * t52198 ) + t51946 * ( t52007 + t51949 * (  - t51955 + t52090 ) + t51947 * t52198 ) ) );
		}
	}
	
protected:
	const Math::Vector< 3, VType >& m_v;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

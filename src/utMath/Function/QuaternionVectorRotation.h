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
 * function class for rotating a 3-vector around a quaternion
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

namespace Ubitrack { namespace Math { namespace Function {
 
/**
 * Function class for rotating a 3-vector around a quaternion
 *
 */
struct QuaternionVectorRotation
{
	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 3; }
	
	/**
	 * Evaluate the function on the input \c input and store the result in
	 * \c result. \c VT1, \c VT2 and \c VT3 can be assumed to be Boost uBlas vectors
	 * using operator() to access elements
	 */
	template< class VT1, class VT2, class VT3 > 
	void evaluate( VT1& result, const VT2& q, const VT3& vec ) const
	{
		// precomputation of some values
		typename VT1::value_type xy = q( 0 ) * q( 1 );
		typename VT1::value_type xz = q( 0 ) * q( 2 );
		typename VT1::value_type yz = q( 1 ) * q( 2 );
		typename VT1::value_type ww = q( 3 ) * q( 3 );
		typename VT1::value_type wx = q( 3 ) * q( 0 );
		typename VT1::value_type wy = q( 3 ) * q( 1 );
		typename VT1::value_type wz = q( 3 ) * q( 2 );

		result( 0 ) = vec( 0 ) * ( 2*(q( 0 )*q( 0 ) + ww) - 1 ) + vec( 1 ) * 2 * (xy - wz) + vec( 2 ) * 2 * (wy + xz);
		result( 1 ) = vec( 0 ) * 2 * (xy + wz) + vec( 1 ) * ( 2*(q( 1 )*q( 1 ) + ww) - 1 ) + vec( 2 ) * 2 * (yz - wx);
		result( 2 ) = vec( 0 ) * 2 * (xz - wy) + vec( 1 ) * 2 * (wx + yz) + vec( 2 ) * ( 2*(q( 2 )*q( 2 ) + ww) - 1 );
	}
	
	/**
	 * Evaluate the function on the input \c input and return both the result
	 * and the jacobian. All parameters can be assumed to be Boost 
	 * uBlas vectors and matrices using operator() to access elements.
	 */
	template< class VT1, class VT2, class VT3, class MT1, class MT2 > 
	void evaluateWithJacobian( VT1& result, const VT2& q, const VT3& vec, MT1& jacQ, MT2& jacVec ) const
	{
		evaluate( result, q, vec );
		jacobian( q, vec, jacQ, jacVec );
	}
	
	/**
	 * Evaluate the jacobian on the input \c input.
	 * All parameters can be assumed to be Boost uBlas vectors and matrices 
	 * using operator() to access elements.
	 */
	template< class VT2, class VT3, class MT1, class MT2 > 
	void jacobian( const VT2& q, const VT3& vec, MT1& jacQ, MT2& jacVec ) const
	{
		if ( jacQ.size1() == 3 )
		{
			typedef typename MT1::value_type VType;
			VType t257190 = q( 1 );
			VType t257188 = vec( 0 );
			VType t257187 = q( 0 );
			VType t257191 = vec( 1 );
			VType t257194 = vec( 2 );
			VType t257193 = q( 2 );
			VType t257199 = q( 3 );
			VType t257189 = 2 * t257187 * t257188;
			VType t257192 = 2 * t257190 * t257191;
			VType t257195 = 2 * t257193 * t257194;
			VType t257196 = t257189 + t257192 + t257195;
			VType t257206 = 2 * t257188 * t257199;
			VType t257207 =  - 2 * t257191 * t257193;
			VType t257208 = 2 * t257190 * t257194;
			VType t257209 = t257206 + t257207 + t257208;
			VType t257214 = 2 * t257188 * t257193;
			VType t257215 = 2 * t257191 * t257199;
			VType t257216 =  - 2 * t257187 * t257194;
			VType t257217 = t257214 + t257215 + t257216;
			VType t257197 =  - 2 * t257188 * t257190;
			VType t257198 = 2 * t257187 * t257191;
			VType t257200 = 2 * t257194 * t257199;
			VType t257201 = t257197 + t257198 + t257200;
			jacQ( 0, 0 ) = t257196;
			jacQ( 0, 1 ) = t257201;
			jacQ( 0, 2 ) =  - 2 * t257188 * t257193 + 2 * t257187 * t257194 - 2 * t257191 * t257199;
			jacQ( 0, 3 ) = t257209;
			jacQ( 1, 0 ) = 2 * t257188 * t257190 - 2 * t257187 * t257191 - 2 * t257194 * t257199;
			jacQ( 1, 1 ) = t257196;
			jacQ( 1, 2 ) = t257209;
			jacQ( 1, 3 ) = t257217;
			jacQ( 2, 0 ) = t257217;
			jacQ( 2, 1 ) = 2 * t257191 * t257193 - 2 * t257190 * t257194 - 2 * t257188 * t257199;
			jacQ( 2, 2 ) = t257196;
			jacQ( 2, 3 ) = t257201;		
		}
		
		if ( jacVec.size1() == 3 )
		{
			typedef typename MT2::value_type VType;
			VType t257336 = q( 1 );
			VType t257339 = q( 2 );
			VType t257345 = q( 3 );
			VType t257343 = q( 0 );
			VType t257340 = t257339 * t257339;
			VType t257341 =  - 2 * t257340;
			VType t257356 = t257343 * t257343;
			VType t257357 =  - 2 * t257356;
			VType t257337 = t257336 * t257336;
			VType t257338 =  - 2 * t257337;
			jacVec( 0, 0 ) = 1 + t257338 + t257341;
			jacVec( 0, 1 ) = 2 * t257336 * t257343 - 2 * t257339 * t257345;
			jacVec( 0, 2 ) = 2 * ( t257339 * t257343 + t257336 * t257345 );
			jacVec( 1, 0 ) = 2 * ( t257336 * t257343 + t257339 * t257345 );
			jacVec( 1, 1 ) = 1 + t257341 + t257357;
			jacVec( 1, 2 ) = 2 * t257336 * t257339 - 2 * t257343 * t257345;
			jacVec( 2, 0 ) = 2 * t257339 * t257343 - 2 * t257336 * t257345;
			jacVec( 2, 1 ) = 2 * ( t257336 * t257339 + t257343 * t257345 );
			jacVec( 2, 2 ) = 1 + t257338 + t257357;		
		}
	}
};

} } } // namespace Ubitrack::Math::Function

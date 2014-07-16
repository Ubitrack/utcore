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
 * Several casts of rotation types 
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */


#ifndef __UBITRACK_MATH_UTIL_ROTATION_CAST_H_INCLUDED__
#define __UBITRACK_MATH_UTIL_ROTATION_CAST_H_INCLUDED__

#include <cmath> // std::sqrt, std::cos, std::sin

namespace Ubitrack { namespace Math { 

//forward declaration to some templated datastructures
class Quaternion;
template< typename T, std::size_t N > class Vector;
template< typename T, std::size_t M, std::size_t N > class Matrix;


namespace Util {

/**
 * @brief Rotation cast functor to cast various rotation representations
 *
 * This functor supports the change of a certain rotation 
 * representation (encoding) from one type into another type.
 *
 * Each specialization can be used to accept all reasonable
 * rotation representations and cast to is own representation.
 *
 * The conversions are mainly inspired from Martin Baker's
 * website: http://www.euclideanspace.com/maths/geometry/rotations/conversions/ 
 *
 * @tparam rotation_type names the desired rotation representation.
 */
 
template< typename rotation_type >
struct RotationCast
{
	typedef rotation_type return_type;
	
	template< typename rotation_in >
	return_type operator()( const rotation_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_ROTATION_CONVERSION_AVAILABLE );
		return rotation_type();
	}
};

/** @internal  specialization to cast a rotation into a quaternion representation. */
template< >
struct RotationCast< Math::Quaternion >
{
	typedef Math::Quaternion return_type;
	
	// function to catch so far non existing rotation casts
	template< typename rotation_in >
	return_type operator()( const rotation_in & rotation ) const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_ROTATION_CONVERSION_AVAILABLE );
		return return_type();
	}

	// if input type equals output type: no change at all
	return_type operator()( const return_type & other )const
	{
		return other;
	}
	
	template< typename T, std::size_t M, std::size_t N >
	return_type operator()( const Matrix< T, M, N > & matrix ) const
	{
		UBITRACK_STATIC_ASSERT( ( ( M == 3 || M == 4 ) && ( N == 3 || N == 4 ) ), NO_MATCHING_MATRIX_REPRESENTATION_FOR_ROTATION );
		const T m00div4 = matrix( 0, 0 ) / 4;
		const T m11div4 = matrix( 1, 1 ) / 4;
		const T m22div4 = matrix( 2, 2 ) / 4;
		const T qw = std::sqrt( std::max( 0., ( 0.25 + m00div4 + m11div4 + m22div4 ) ) ) ; 
		const T qxt = std::sqrt( std::max( 0., ( 0.25 + m00div4 - m11div4 - m22div4 ) ) ) ; 
		const T qyt = std::sqrt( std::max( 0., ( 0.25 - m00div4 + m11div4 - m22div4 ) ) ) ; 
		const T qzt = std::sqrt( std::max( 0., ( 0.25 - m00div4 - m11div4 + m22div4 ) ) ) ; 
		const T qx = ( matrix( 2, 1 ) - matrix( 1, 2 ) ) < 0 ? -qxt : qxt;
		const T qy = ( matrix( 0, 2 ) - matrix( 2, 0 ) ) < 0 ? -qyt : qyt;
		const T qz = ( matrix( 1, 0 ) - matrix( 0, 1 ) ) < 0 ? -qzt : qzt;
		return return_type( qx, qy, qz, qw );
	}
	
	template< typename T >
	return_type operator()( const Math::Vector< T, 4 > & axisAngle ) const
	{
		const T angle = axisAngle[ 3 ];
		const T halfAngle = angle / 2;
		const T sinHalfAngle = std::sin( halfAngle );
		
		const T qx = axisAngle[ 0 ] * sinHalfAngle;
		const T qy = axisAngle[ 1 ] * sinHalfAngle;
		const T qz = axisAngle[ 2 ] * sinHalfAngle;
		const T qw = std::cos( halfAngle );
		return return_type( qx, qy, qz, qw );
	}
	
	template< typename T >
	return_type operator()( const Math::Vector< T, 3 > & axisAngle ) const
	{
		const T x = axisAngle[ 0 ];
		const T y = axisAngle[ 1 ];
		const T z = axisAngle[ 2 ];
		const T angle = std::sqrt( x*x+y*y+z*z );
		if ( angle < 1e-10 )
			return return_type( 0, 0, 0, 1 );

		const T rx = axisAngle[ 0 ] / angle;
		const T ry = axisAngle[ 1 ] / angle;
		const T rz = axisAngle[ 2 ] / angle;
		const Math::Vector< T, 4 > aa( rx, ry, rz, angle );
		return RotationCast< return_type >( )( aa );
	}

};

/** @internal specialization to cast any rotation into a axis/angle( rx, ry, rz, angle ) representation. */
template< typename T >
struct RotationCast< Math::Vector< T, 4 > >
{
	typedef Math::Vector< T, 4 > return_type;
	
	// function to catch so far non existing rotation casts
	template< typename rotation_in >
	return_type operator()( const rotation_in & rotation )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_ROTATION_CONVERSION_AVAILABLE );
		return return_type();
	}

	// if input type equals output type: no change at all
	return_type operator()( const return_type & other ) const
	{
		return other;
	}
	
	return_type operator()( const Math::Quaternion & quat ) const
	{
		const T angle = 2 * std::acos( quat.w() );
		const T divisor = std::sqrt( 1 - quat.w()*quat.w() );
		const T x = quat.x() / divisor;
		const T y = quat.y() / divisor;
		const T z = quat.z() / divisor;
		return Math::Vector< T, 4 >( x, y, z, angle );
	}
	
	return_type operator()( const Math::Vector< T, 3 > & rotAxis ) const
	{
		const T x = rotAxis[ 0 ];
		const T y = rotAxis[ 1 ];
		const T z = rotAxis[ 2 ];
		const T angle = std::sqrt( x*x+y*y+z*z );
		if( angle != angle || angle < 1e-10 )
			return return_type( 0, 0, 0, 0 );
			
		return return_type( x/angle, y/angle, z/angle, angle );
	}
	
	
	// commented out since it does not work as is: see tests please
	template< typename ValType, std::size_t M, std::size_t N >
	return_type operator()( const Matrix< ValType, M, N > & matrix ) const
	{
		UBITRACK_STATIC_ASSERT( ( ( M == 3 || M == 4 ) && ( N == 3 || N == 4 ) ), NO_MATCHING_MATRIX_REPRESENTATION_FOR_ROTATION );
		
		//easy way here
		// const Math::Quaternion quat = RotationCast< Math::Quaternion >()( matrix );
		// return RotationCast< return_type >()( quat );
		
		/// @todo put in more checks, this transformation is not really reliable at the moment
		const T angle = std::acos( ( matrix( 0, 0 ) + matrix( 1, 1 ) + matrix( 2, 2 ) - 1 ) / 2 );
		if ( angle != angle || angle < 1e-10 )
			return return_type( 0, 0, 0, 0 );
		
		//std::cout << "Angle " << angle ;
		const T m2112 = matrix( 2, 1 ) - matrix( 1, 2 );
		const T m0220 = matrix( 0, 2 ) - matrix( 2, 0 );
		const T m1001 = matrix( 1, 0 ) - matrix( 0, 1 );
		
		const T divisor = sqrt( m2112*m2112 + m0220*m0220 + m1001*m1001 );
		//std::cout << " Divisor " << divisor << "\n";
		const T rx = m2112 / divisor;
		const T ry = m0220 / divisor;
		const T rz = m1001 / divisor;
		return return_type( rx, ry, rz, angle );
		
		// alternative way:
		// const T divisor = 2*std::sin( angle ) ;//sqrt( m2112*m2112 + m0220*m0220 + m1001*m1001 );
		// //std::cout << " Divisor " << divisor << "\n";
		// const T rx = m2112 / divisor;
		// const T ry = m0220 / divisor;
		// const T rz = m1001 / divisor;
		// return return_type( rx, ry, rz, angle );
	}
};

/**
 * @internal specialization to cast any rotation into a axis/angle(rx, ry, rz)
 * representation where the norm encodes the rotation angle.
 */
template< typename T >
struct RotationCast< Math::Vector< T, 3 > >
{
	typedef Math::Vector< T, 3 > return_type;
	
	template< typename rotation_in >
	return_type operator()( const rotation_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_ROTATION_CONVERSION_AVAILABLE );
		return return_type();
	}
	
	// if input type equals output type: no change at all
	return_type operator()( const return_type & other )const
	{
		return other;
	}
	
	return_type operator()( const Math::Vector< T, 4 > & axisAngle ) const
	{
		const T x = axisAngle[ 0 ];
		const T y = axisAngle[ 1 ];
		const T z = axisAngle[ 2 ];
		const T a = axisAngle[ 3 ];
		const T dst = std::sqrt( x*x+y*y+z*z );
		if( dst != dst || dst < 1e-10 )
			return Math::Vector< T, 3 >( 0, 0, 0 );
			
		return Math::Vector< T, 3 >( a*x/dst, a*y/dst, a*z/dst );
	}
	
	return_type operator()( const Math::Quaternion & quat ) const
	{
		const Math::Vector< T, 4 > axisAngle = RotationCast< Math::Vector< T, 4 > >()( quat );
		return RotationCast< return_type >()( axisAngle );
	}
	
	template< typename ValType, std::size_t M, std::size_t N >
	return_type operator()( const Math::Matrix< ValType, M, N > & matrix ) const
	{
		const Math::Vector< T, 4 > axisAngle = RotationCast< Math::Vector< T, 4 > >()( matrix );
		return RotationCast< return_type >()( axisAngle );
	}
};

/** @internal specialization to cast any rotation into a matrix representation. */
template< typename T, std::size_t M, std::size_t N >
struct RotationCast< Matrix< T, M, N > >
{
	typedef Matrix< T, M, N > return_type;
	
	RotationCast()
	{
		UBITRACK_STATIC_ASSERT( ( ( M == 3 || M == 4 ) && ( N == 3 || N == 4  ) ), NO_MATCHING_MATRIX_REPRESENTATION_FOR_ROTATION );
	};
	
	// function to catch so far non existing rotation casts
	template< typename rotation_in >
	return_type operator()( const rotation_in & rotation )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_ROTATION_CONVERSION_AVAILABLE );
		return return_type();
	}

	// if input type equals output type: no change at all
	return_type operator()( const return_type & other )const
	{
		return other;
	}
	
	return_type operator()( const Math::Quaternion & other ) const
	{
		const T ww = other.w() * other.w();
		const T xx = other.x() * other.x();
		const T yy = other.y() * other.y();
		const T zz = other.z() * other.z();
		const T xy = other.x() * other.y();
		const T yz = other.y() * other.z();
		const T zx = other.z() * other.x();
		const T wx = other.w() * other.x();
		const T wy = other.w() * other.y();
		const T wz = other.w() * other.z();
		
		return_type result;
		// first matrix row
		result( 0, 0 ) = ww + xx - yy - zz;
		result( 1, 0 ) = 2* ( xy + wz );
		result( 2, 0 ) = 2* ( zx - wy );
		// second matrix row
		result( 0, 1 ) = 2* ( xy - wz );
		result( 1, 1 ) =  ww - xx + yy - zz;
		result( 2, 1 ) = 2* ( yz + wx );
		// third matrix row
		result( 0, 2 ) = 2* ( zx + wy );
		result( 1, 2 ) = 2* ( yz - wx );
		result( 2, 2 ) = ww - xx - yy + zz;
		return result;
			
		// T values [ 9 ];
		// //following row_major order:
		// values[ 0 ] = ww + xx - yy - zz;
		// values[ 1 ] = 2* ( xy - wz );
		// values[ 2 ] = 2* ( zx + wy );
		// values[ 3 ] = 2* ( xy + wz );
		// values[ 4 ] =  ww - xx + yy - zz;
		// values[ 5 ] = 2* ( yz - wx );
		// values[ 6 ] = 2* ( zx - wy );
		// values[ 7 ] = 2* ( yz + wx );
		// values[ 8 ] =  ww - xx - yy + zz;
		// return return_type	( values );
	}
	
	
	return_type operator()( const Math::Vector< T, 4 > & axisAngle ) const
	{
		const T rx = axisAngle[ 0 ];
		const T ry = axisAngle[ 1 ];
		const T rz = axisAngle[ 2 ];

		const T c = std::cos( axisAngle[ 3 ] );
		const T s = std::sin( axisAngle[ 3 ] );
		const T t = 1 - c;

		return_type matrix;
		// diagonal entries first
		matrix( 0, 0 ) = c + rx*rx*t;
		matrix( 1, 1 ) = c + ry*ry*t;
		matrix( 2, 2 ) = c + rz*rz*t;

		const T tmp1 = rx*ry*t;
		const T tmp2 = rz*s;
		matrix( 1, 0 ) = tmp1 + tmp2;
		matrix( 0, 1 ) = tmp1 - tmp2;
		
		const T tmp3 = rz*rx*t;
		const T tmp4 = ry*s;
		matrix( 2, 0 ) = tmp3 - tmp4;
		matrix( 0, 2 ) = tmp3 + tmp4;
		
		const T tmp5 = ry*rz*t;
		const T tmp6 = rx*s;
		matrix( 2, 1 ) = tmp5 + tmp6;
		matrix( 1, 2 ) = tmp5 - tmp6;
		return matrix;
	}
	
	return_type operator()( const Math::Vector< T, 3 > &rotAxis ) const
	{
		const Math::Vector< T, 4 > axisAngle = RotationCast< Math::Vector< T, 4 > >()( rotAxis );
		return RotationCast< return_type >() ( axisAngle );
	}
};

} } } //namespace Ubitrack::Math::Util

#endif // __UBITRACK_MATH_UTIL_ROTATION_CAST_H_INCLUDED__
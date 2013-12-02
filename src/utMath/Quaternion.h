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
 * Wrapper for boost::math::quaternion.
 * @author Florian Echtler <echtler@in.tum.de>
 */


#ifndef _Ubitrack_Math_Quaternion_INCLUDED_
#define _Ubitrack_Math_Quaternion_INCLUDED_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/quaternion.hpp>

// WARNING: all boost/serialization headers should be
//          included AFTER all boost/archive headers
#include <boost/serialization/access.hpp>

#include <utCore.h>
#include "Vector.h"
#include <math.h>

namespace Ubitrack { namespace Math {

//forward declaration to matrix template
template< std::size_t M, std::size_t N, typename T > class Matrix;

/**
 * @ingroup math
 * Wraps a boost::math::quaternion for convenience.
 * Provides axis-angle and matrix constructors as well 
 * as some useful operators, e.g. normalization.
 */

class UBITRACK_EXPORT Quaternion
 	: public boost::math::quaternion<double>
{
	friend class boost::serialization::access;

	public:

		/**
		 * default constructor, sometimes needed
		 */
		Quaternion()
			: boost::math::quaternion< double >( 1, 0, 0, 0 )
		{}
		
		Quaternion( double x, double y, double z)			
		{
			double c1 = cos(y/2);
			double s1 = sin(y/2);
			double c2 = cos(x/2);
			double s2 = sin(x/2);
			double c3 = cos(z/2);
			double s3 = sin(z/2);
			double c1c2 = c1*c2;
			double s1s2 = s1*s2;
			a =c1c2*c3 - s1s2*s3;
			b =c1c2*s3 + s1s2*c3;
			c =s1*c2*c3 + c1*s2*s3;
			d =c1*s2*c3 - s1*c2*s3;
		}

		/**
		 * axis-angle constructor
		 * @param axis rotation axis
		 * @param angle rotation angle
		 */
		Quaternion(const Vector<3>& axis, const double angle)
		{
			double length = ::boost::numeric::ublas::norm_2(axis);
			Vector<3> normAxis = axis*(sin(angle/2.0)/length);

			a = cos(angle/2.0);
			b = normAxis(0);
			c = normAxis(1);
			d = normAxis(2);
		}

		/**
		 * value constructor
		 * @param _x first imaginary part
		 * @param _y second imaginary part
		 * @param _z third imaginary part
		 * @param _w real part
		 */
		Quaternion( double _x, double _y, double _z, double _w)
			: boost::math::quaternion< double >( _w, _x, _y, _z )
		{ }

		/**
		 * base class constructor
		 * @param q a boost::math::quaternion to copy
		 */
		Quaternion( const ::boost::math::quaternion< double >& q )
			: boost::math::quaternion< double >( q )
		{ }

		/**
		 * rotation matrix constructor
		 * @param mat a rotation matrix
		 */
		Quaternion( const Math::Matrix< 0, 0, double >& mat );

		/**
		 * get the real part
		 * @return real component w (member a in the base class)
		 */
		double w() const { return a; }
		
		double x() const { return b; }
		double y() const { return c; }
		double z() const { return d; }

		/**
		 * normalize quaternion to describe a rotation
		 * @return the newly normalized quaternion
		 */
		Quaternion& normalize() {
			*this /= boost::math::abs(*this);
			return *this;
		}

		/**
		 * invert quaternion
		 * @result the newly inverted quaternion
		 */
		Quaternion& invert() {
			*this = boost::math::conj(*this);
			return *this;
		}

		/**
		 * return inverted copy of quaternion
		 * @return inverted copy
		 */
		Quaternion operator~() const {
			return Quaternion( boost::math::conj( *this ) );
		}

		/**
		 * rotate a vector by a quaternion
		 * @param x input vector
		 * @return output vector
		 */
		Vector< 3 > operator*( const Vector< 3 >& x ) const;

		/**
		 * @return the angle of the rotation represented by the quaternion
		 */
		double angle() const
		{ return 2 * acos( w() > 0 ? w() : -w() ); }

		/** 
		 * sets the rotation of a given matrix to the value of the quaternion
		 * the rotation will be stored in the upper 3 x 3 subrange of the matrix
		 * @param matrix a boost::numeric::ublas::matrix, where the rotation of the quaternion will be set
		 */
		template< class M >
		void toMatrix( M& matrix ) const;

		/**
		 * sets the given axis-angle parameters to represent the value of this quaternion.
		 * Note: Call normalize() in case of doubt for proper operation.
		 */
		template< class M >
		void toAxisAngle( Math::Vector< 3, M >& axis, M& angle );
		
		/**
		 * computes the quaternion logarithm
		 */
		Vector< 3 > toLogarithm() const;
		
		/**
		 * creates a quaternion from a quaternion logarithm
		 */
		static Quaternion fromLogarithm( const Vector< 3 >& v );

		/**
		 * negates the quaternion if the result is closer to a given reference.
		 * Necessary for analysis of quaternion sequences.
		 * @param ref reference quaternion
		 * @return original or negated quaternion
		 */
		Quaternion negateIfCloser( const Quaternion& ref ) const
		{
			double prod( x() * ref.x() + y() * ref.y() + z() * ref.z() + w() * ref.w() );
			if ( prod >= 0 )
				return *this;
			else
				return Quaternion( -x(), -y(), -z(), -w() );
		}

		/**
		 * convert quaternion to Euler angles in the common z-y-x order (global reference system)
		 * @return rx,ry,rz result angles
		 * ### deprecated
		 */
		Vector<3> getEulerAngles() const;

		/**
		 * Sequence of rotations for computation of Euler angles. XYZ
		 * means rotation about X comes before rotation about Y before rotation about Z.
		 */
		typedef enum {
		     EULER_SEQUENCE_XYZ = 0,
			 EULER_SEQUENCE_YZX,
			 EULER_SEQUENCE_ZXY,
			 EULER_SEQUENCE_ZYX,
			 EULER_SEQUENCE_XZY,
			 EULER_SEQUENCE_YXZ
		} t_EulerSequence;

		/**
		 * Convert quaternion to Euler angles in the given sequence.
		 * @return rx,ry,rz result angles
		 */
		Vector<3> getEulerAngles( t_EulerSequence seq ) const;
		

		/**
		 * Store the quaternion in a Boost uBlas vector.
		 * The order is x, y, z, w (as usual).
		 */
		template< class VType >
		void toVector( boost::numeric::ublas::vector_expression< VType >& v ) const
		{ 
			typedef typename VType::value_type T; 
			v()( 0 ) = T( x() ); 
			v()( 1 ) = T( y() ); 
			v()( 2 ) = T( z() ); 
			v()( 3 ) = T( w() ); 
		}
		
		/**
		 * Retrieve a quaternion from a Boost uBlas vector.
		 * The order is x, y, z, w (as usual).
		 */
		template< class VType >
		static Quaternion fromVector( const VType& v )
		{ return Quaternion( v( 0 ), v( 1 ), v( 2 ), v( 3 ) ); }
		
	protected:

		/// serialize Quaternion object
		template< class Archive > 
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & b;
			ar & c;
			ar & d;
			ar & a;
		}

};


/// stream output operator
UBITRACK_EXPORT std::ostream& operator<<( std::ostream& s, const Quaternion& q );

template< class M >
void Quaternion::toMatrix( M& matrix ) const
{
	typedef typename M::value_type T;

	T X = T( x() );
	T Y = T( y() );
	T Z = T( z() );
	T W = T( w() );

	T xx = X * X;
	T xy = X * Y;
	T xz = X * Z;
	T xw = X * W;
	T yy = Y * Y;
	T yz = Y * Z;
	T yw = Y * W;
	T zz = Z * Z;
	T zw = Z * W;

	matrix(0,0) = 1 - 2 * ( yy + zz );
	matrix(0,1) =     2 * ( xy - zw );
	matrix(0,2) =     2 * ( xz + yw );
	matrix(1,0) =     2 * ( xy + zw );
	matrix(1,1) = 1 - 2 * ( xx + zz );
	matrix(1,2) =     2 * ( yz - xw );
	matrix(2,0) =     2 * ( xz - yw );
	matrix(2,1) =     2 * ( yz + xw );
	matrix(2,2) = 1 - 2 * ( xx + yy );
}


// Code taken from http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm
template< class M >
void Quaternion::toAxisAngle( Math::Vector< 3, M >& axis, M& angle )
{
	// assuming quaternion normalised then w is less than 1, so term always positive.
	angle = 2 * acos( a );
	
	double s = sqrt( 1 - a*a );
	// test to avoid divide by zero, s is always positive due to sqrt
	if (s < 0.001) { 
		// if s close to zero then direction of axis not important
		axis( 0 ) = b; // if it is important that axis is normalised then replace with x=1; y=z=0;
		axis( 1 ) = c;
		axis( 2 ) = d;
	} else {
		axis( 0 ) = b / s; // normalise axis
		axis( 1 ) = c / s;
		axis( 2 ) = d / s;
	}
}


/**
 * SLERP
 * spherical linear interpolation between two quaternions
 *
 * @param a first quaternion
 * @param b second quaternion
 * @param t interpolation point, 0.0 represents a, 1.0 represents b
 * @return result quaternion
 */
UBITRACK_EXPORT Quaternion slerp( const Quaternion& a, const Quaternion& b, const double t );

/** same as slerp */
UBITRACK_EXPORT Quaternion linearInterpolate( const Quaternion& a, const Quaternion& b, const double t );


} } // namespace Ubitrack::Math

#endif // _Ubitrack_Math_Quaternion_INCLUDED_


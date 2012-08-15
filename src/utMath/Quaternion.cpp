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


#include "Quaternion.h"
#include <iomanip>
#include "Matrix.h"
#include <boost/math/constants/constants.hpp>


namespace Ubitrack { namespace Math {


// create quaternion from (assumed) rotation matrix
Quaternion::Quaternion( const boost::numeric::ublas::matrix<double>& mat )
{
	double S, X, Y, Z, W;
	double T = 1.0 + mat(0,0) + mat(1,1) + mat(2,2);

	if (T > 0) {
		S = sqrt(T) * 2;
		X = ( mat(1,2) - mat(2,1) ) / S;
		Y = ( mat(2,0) - mat(0,2) ) / S;
		Z = ( mat(0,1) - mat(1,0) ) / S;
		W = 0.25 * S;
	} else if ( (mat(0,0) > mat(1,1)) && (mat(0,0) > mat(2,2)) ) { // Column 0
		S = sqrt( 1.0 + mat(0,0) - mat(1,1) - mat(2,2) ) * 2;
		X = 0.25 * S;
		Y = ( mat(0,1) + mat(1,0) ) / S;
		Z = ( mat(2,0) + mat(0,2) ) / S;
		W = ( mat(1,2) - mat(2,1) ) / S;
	} else if ( mat(1,1) > mat(2,2) ) { // Column 1
		S = sqrt( 1.0 + mat(1,1) - mat(0,0) - mat(2,2) ) * 2;
		X = ( mat(0,1) + mat(1,0) ) / S;
		Y = 0.25 * S;
		Z = ( mat(1,2) + mat(2,1) ) / S;
		W = ( mat(2,0) - mat(0,2) ) / S;
	} else { // Column 2
		S = sqrt( 1.0 + mat(2,2) - mat(0,0) - mat(1,1) ) * 2;
		X = ( mat(2,0) + mat(0,2) ) / S;
		Y = ( mat(1,2) + mat(2,1) ) / S;
		Z = 0.25 * S;
		W = ( mat(0,1) - mat(1,0) ) / S;
	}

	// Quaternion has to be inverted, as the above code from the
	// Matrix and Quaternion FAQ expects a column-major matrix.
	a = -W;
	b =  X;
	c =  Y;
	d =  Z;
	
	normalize();
}


// retrieve Euler angles corresponding to whatever sequence the bloody Kuka robot uses
//
// return value: Vector<3>( angle_x_axis, angle_y_axis, angle_z_axis )
//
// note that the angles are returned in x-y-z order, but are probably applied as z-y-x
// (though nobody can be really sure - ask three times, you'll get four different answers,
// none of them correct). this sequence has been arrived at through a ton of guesswork,
// but it gives proper results when applied to the robot.
Vector<3> Quaternion::getEulerAngles() const {

	double rx = 0.0;
	double ry = 0.0;
	double rz = 0.0;

	//boost::numeric::ublas::matrix<double> m(3,3);
	Matrix<3,3,double> m;
	toMatrix( m );

	//std::cout << m << std::endl;
	//std::cout << std::setprecision(20);
	//std::cout << "m02: " << m(0,2) << std::endl;
	
	double sy = -m(2,0);              // sine value sanity check
	if (sy >  1.0) sy =  1.0;
	if (sy < -1.0) sy = -1.0;

	ry = asin( sy );                  // get y-axis angle
	double cy = cos( ry );
	double tr_x, tr_y;

	if ( fabs( cy ) > 10e-6 ) {       // singularity (ry == pi/2)?

		tr_x =  m(2,2) / cy;            // get x-axis angle
		tr_y =  m(2,1) / cy;

		rx   = atan2( tr_y, tr_x );

		tr_x =  m(0,0) / cy;            // get z-axis angle
		tr_y =  m(1,0) / cy;

		rz   = atan2( tr_y, tr_x );

	} else {                          // singularity (aka gimbal lock)

		rx   = 0;                       // assume x-axis angle as zero

		tr_x =  m(1,1);                 // get z-axis angle
		tr_y = -m(0,1);

		rz   = atan2( tr_y, tr_x );
	}

	return Vector<3>( rx, ry, rz ); 
}


/*
 * This implementation is based on:
 * http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/index.htm
 * Note in particular the section about "Alternative Euler Angle
 * Sequences" and the referenced literature (available for download).
 */
Vector<3> Quaternion::getEulerAngles( Quaternion::t_EulerSequence seq ) const
{
	int sign;

	double x_,y_,z_;
	double w_ = w();

	if ( seq == EULER_SEQUENCE_XYZ ) 
	{
		x_ = x(); y_ = y(); z_ = z();
		sign = -1;
	}
	else if ( seq == EULER_SEQUENCE_YZX ) 
	{
		x_ = y(); y_ = z(); z_ = x();
		sign = -1;
	}
	else if ( seq == EULER_SEQUENCE_ZXY )
	{
		x_ = z(); y_ = x(); z_ = y();
		sign = -1;
	}
	else if ( seq == EULER_SEQUENCE_ZYX )
	{
		x_ = z(); y_ = y(); z_ = x();
		sign =  1;
	}
	else if ( seq == EULER_SEQUENCE_XZY )
	{
		x_ = x(); y_ = z(); z_ = y();
		sign =  1;
	}
	else // ( seq == EULER_SEQUENCE_YXZ )
	{
		x_ = y(); y_ = x(); z_ = z();
		sign =  1;
	}

	double alpha, gamma;
	double beta  = 2 * ( w_*y_ + sign*x_*z_ );

	// Singularity at north pole
	if ( beta > 0.998 ) 
	{
		alpha = 2 * atan2( x_, w_ );
		beta  = boost::math::constants::pi<double>() / 2;
		gamma = 0;

		return Vector<3>( alpha, beta, gamma );
	}
	// Singularity at south pole
	if ( beta < -0.998 ) 
	{
		alpha = -2 * atan2( x_, w_ );
		beta  = -boost::math::constants::pi<double>() / 2;
		gamma = 0;

		return Vector<3>( alpha, beta, gamma );
	}

	alpha = atan2( 2*( w_*x_ - sign*y_*z_ ), ( 1 - 2*(x_*x_ + y_*y_ ) ) );
	beta  = asin( beta );
	gamma = atan2( 2*( w_*z_ - sign*x_*y_ ), ( 1 - 2*(y_*y_ + z_*z_ ) ) );

	return Vector<3>( alpha, beta, gamma ); 
}


// rotate a vector by a quaternion
Vector< 3 > Quaternion::operator*( const Vector< 3 >& vec ) const
{
	Vector< 3 > r;

	// precomputation of some values
	double xy = x() * y();
	double xz = x() * z();
	double yz = y() * z();
	double ww = w() * w();
	double wx = w() * x();
	double wy = w() * y();
	double wz = w() * z();

	r( 0 ) = vec( 0 ) * ( 2*(x()*x() + ww) - 1 ) + vec( 1 ) * 2 * (xy - wz) + vec( 2 ) * 2 * (wy + xz);
	r( 1 ) = vec( 0 ) * 2 * (xy + wz) + vec( 1 ) * ( 2*(y()*y() + ww) - 1 ) + vec( 2 ) * 2 * (yz - wx);
	r( 2 ) = vec( 0 ) * 2 * (xz - wy) + vec( 1 ) * 2 * (wx + yz) + vec( 2 ) * ( 2*(z()*z() + ww) - 1 );

	return r;

}


Vector< 3 > Quaternion::toLogarithm() const
{
	// always take the quaternion with w > 0
	double s = w() >= 0 ? 1 : -1;
	double omega;
	if ( w() * s < 1.0 )
		omega = 2 * acos( w() * s );
	else
		omega = 0.0;
	
	double imagLen = sqrt( x()*x() + y()*y() + z()*z() );
	if ( imagLen > 1e-12 )
	{
		s = s * omega / imagLen;
		return Vector< 3 >( x() * s, y() * s, z() * s );
	}
	else
		return Vector< 3 >( 0, 0, 0 );
}


Quaternion Quaternion::fromLogarithm( const Vector< 3 >& v )
{
	double omega = boost::numeric::ublas::norm_2( v );
	if ( omega > 1e-12 )
	{
		double s = sin( omega / 2 ) / omega;
		return Quaternion( s * v( 0 ), s * v( 1 ), s * v( 2 ), cos( omega / 2 ) );
	}
	else
		return Quaternion( 0, 0, 0, 1 );
}


// stream output operator
std::ostream& operator<<( std::ostream& s, const Quaternion& q )
{
	s << "[ ( " << q.x() << " "
	            << q.y() << " "
	            << q.z() << " ) "
	            << q.w() << " ]";
	return s;
}


// spherical linear interpolation between two quaternions
Quaternion slerp (const Quaternion& a, const Quaternion& b, const double t)
{
  // make working copies of the input parameters
  Quaternion x = a;
  Quaternion y = b;

  // calculate the norm of the difference between a and b
  Quaternion difference = x-y;
  double diffNorm = boost::math::norm(difference);

  // if norm is too large, we negate the first quaternion
  // note that this does not change the rotation which is represented;
  if (diffNorm > 2.0)
  {
    x = Quaternion ( -a.x(), -a.y(), -a.z(), -a.w() );
  }

  // calculate the dot product
  double dotProduct = x.x()*y.x()+x.y()*y.y()+x.z()*y.z()+x.w()*y.w();

  // calculate weights w1 and w2 for a linear combination of x and y;
  double w1, w2;

  // if the angle between x and y is too small we perform a simple linear
  // interpolation
  if ( dotProduct > 0.9999 )
  {
    // small angle => linear interpolation
    w1 = 1.0 - t;
    w2 = t;
  }
  else
  {
    // calucate omega, sin(omega)
    double omega = acos ( dotProduct );
    double sinOmega = sin ( omega );

    // calculate the weights w1 and w2
    w1 = sin( ( 1.0 - t ) * omega ) / sinOmega;
    w2 = sin( t * omega ) / sinOmega;
  }

  // calculate result as weighted sum of x and y;
  Quaternion result( w1*x + w2*y );

  // normalize again for good measure;
  // probably unnecessary
  result.normalize();

  return result;
}


Quaternion linearInterpolate( const Quaternion& a, const Quaternion& b, const double t )
{
	return slerp( a, b, t );
}


} } // namespace Ubitrack::Math


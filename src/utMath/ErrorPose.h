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
 * 6D pose with associated 6x6 covariance matrix
 * @author Daniel Pustka <pustka@in.tum.de>
 */


#ifndef __UBITRACK_MATH_ERRORPOSE_H_INCLUDED__
#define __UBITRACK_MATH_ERRORPOSE_H_INCLUDED__

#include <utCore.h>
// #include "Pose.h"
// #include "Matrix.h"
#include "ErrorVector.h"

namespace Ubitrack { namespace Math {

//forward declaration to reduce header including
class Pose;
template< typename T, std::size_t M, std::size_t N > class Matrix;

/**
 * @ingroup math
 * Represents a 6D pose with errors, composed of a quaternion, a 3-dimensional translation vector
 * and a 6x6 covariance matrix.
 *
 * A pose with errors is a transformation of vectors x in 3D, defined by
@verbatim
x' = q * e_r * x * e_r^* * q^* + t + e_t
@endverbatim
 * where e_r is the error in orientation, represented as a small quaternion e_r=(e_rx, e_ry, e_rz, 1) 
 * and e_t=(e_tx, e_ty, e_tz) is the error in translation.
 * 
 * The covariance matrix of a pose describes the distribution of the 6-vector
 * (e_tx, e_ty, e_tz, e_rx, e_ry, e_rz) around 0.
 */
class UBITRACK_EXPORT ErrorPose
	: public Pose
{

	protected:

	/** the 6-by-6 matrix that represents the covariance */
	Matrix< double, 6, 6 > m_covariance;
	
	public:

		typedef double value_type;

		/** doesn't make much sense, but sometimes we need a default constructor */
		ErrorPose()
		{}

		/**
		 * construct from quaternion, translation and error
		 * @param q a rotation
		 * @param t a translation
		 * @param c a 6x6 covariance matrix
		 */
		template< class MT >
		ErrorPose( const Quaternion& q, const Vector< double, 3 >& t, const MT& c )
			: Pose( q, t )
			, m_covariance( c )
		{}

		/**
		 * construct from pose and error
		 * @param p a pose
		 * @param c a 6x6 covariance matrix
		 */
		template< class MT >
		ErrorPose( const Pose& p, const MT& c )
			: Pose( p )
			, m_covariance( c )
		{}

		/**
		 * get the covariance
		 * @return covariance of the pose
		 */
		const Matrix< double, 6, 6 >& covariance() const
		{ return m_covariance; }

		/**
		 * inverts a pose with errors
		 * @return the inverted error pose
		 */
		ErrorPose operator~( ) const;

		/**
		 * converts the pose to a ublas vector and the multiplicative 6x6 covariance to an 
		 * additive 7x7 covariance matrix.
		 */
		void toAdditiveErrorVector( ErrorVector< double, 7 >& v );

		/**
		 * creates an ErrorPose from a ublas vector with an additive 7x7 covariance.
		 */
		static ErrorPose fromAdditiveErrorVector( const ErrorVector< double, 7 >& v );

	protected:

		friend class ::boost::serialization::access;
		
		/// serialize Pose object
		template< class Archive > 
		void serialize( Archive& ar, const unsigned int version )
		{
			Pose::serialize( ar, version );
			ar & m_covariance;
		}
};


/**
 * @internal multiplies two poses and propagates the error of both poses
 */
UBITRACK_EXPORT ErrorPose operator*( const ErrorPose& a, const ErrorPose& b );

/**
 * @internal multiplies two poses and propagates the error of the first pose
 */
UBITRACK_EXPORT ErrorPose operator*( const ErrorPose& a, const Pose& b );

/**
 * @internal multiplies two poses and propagates the error of the second pose
 */
UBITRACK_EXPORT ErrorPose operator*( const Pose& a, const ErrorPose& b );

/**
 * @internal 
 * Propagates the error to a point of interest.
 * The function performs a pose-times-vector multiplication from a coordinate frame X to a coordinate frame Y.
 * @param a pose with error
 * @param b point in the coordinate frame X described by \ a
 * @return point in coordinate frame Y with associated 3x3-covariance
 */
UBITRACK_EXPORT ErrorVector< double, 3 > operator*( const ErrorPose& a, const Math::Vector< double, 3 >& b );

/**
 * Multiplies two poses and inverts one ( A^-1 * B ).
 * Gives more realistic error propagation if computed together.
 */
UBITRACK_EXPORT ErrorPose invertMultiply( const ErrorPose& a, const ErrorPose& b );

UBITRACK_EXPORT ErrorVector< double, 3 > operator*( const ErrorPose& a, const Math::ErrorVector< double, 3 >& b );

/**
 * performs a linear interpolation between two poses
 * using SLERP and vector interpolation
 * @param x first pose
 * @param y second pose
 * @param t interpolation point between 0.0 and 1.0
 * @return the new pose
 */
UBITRACK_EXPORT ErrorPose linearInterpolate ( const ErrorPose& x, const ErrorPose& y, double t );

/// @internal stream output operator
UBITRACK_EXPORT std::ostream& operator<<( std::ostream& s, const ErrorPose& p );

} } // namespace Ubitrack::Math

#endif // __UBITRACK_MATH_ERRORPOSE_H_INCLUDED__


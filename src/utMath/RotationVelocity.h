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
 * defines a rotation of constant velocity
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_MATH_ROTATIONVELOCITY_H_INCLUDED__
#define __UBITRACK_MATH_ROTATIONVELOCITY_H_INCLUDED__

#include "Vector.h"
#include "Quaternion.h"

namespace Ubitrack { namespace Math {

/**
 * Describes rotation of constant velocity around a single axis.
 *
 * Contains a 3-vector where the norm corresponds to the angular velocity in
 * radians/s, and the direction of the vector describes the rotation axis.
 */
class RotationVelocity
	: public Math::Vector< 3 >
{
public:
	/** default constructor */
	RotationVelocity()
	{}

	/**
	 * Construct from vector_expression
	 * @param e a vector_expression
	 */
	template< class AE > 
	RotationVelocity( const boost::numeric::ublas::vector_expression< AE >& e )
		: Math::Vector< 3 >( e )
	{}
	
	/**
	 * constructs the velocity from three rotation speeds around the x, y and z axes
	 */
	RotationVelocity( double x, double y, double z )
		: Math::Vector< 3 >( x, y, z )
	{}
	
	/**
	 * Construct the velocity from two measurements and the time difference.
	 * Assumes a motion of @verbatim b = a * delta^interval @endverbatim and computes delta.
	 */
	RotationVelocity( const Quaternion& a, const Quaternion& b, double interval )
	{
		Quaternion delta = ~ ( a.negateIfCloser( b ) ) * b;
		*this = delta.toLogarithm() / interval;
	}

	/**
	 * Integrates the velocity over the given time and returns the resulting rotation quaternion.
	 * @param interval integration time in seconds
	 */
	Quaternion integrate( double interval ) const
	{
		return Quaternion::fromLogarithm( Math::Vector< 3 >( *this * interval ) );
	}

	/**
	 * return the total angular velocity in radians/s
	 */
	double angularVelocity() const
	{
		return boost::numeric::ublas::norm_2( *this );
	}
	
	/**
	 * the axis of rotation (only valid if \c totalVelocity > 0)
	 * @return the normalized rotation axis.
	 */
	Math::Vector< 3 > axis() const
	{
		double norm = angularVelocity();
		if ( norm != 0 )
			return *this / norm;
		else
			return *this;
	}
	
};

} } // namespace Ubitrack::Math

#endif

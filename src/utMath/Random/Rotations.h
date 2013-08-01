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
 * Implements functors for generating randomly distributed rotations.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __H__RANDOM_ROTATIONS_H__
#define __H__RANDOM_ROTATIONS_H__
//std
#include <functional>

//boost
#include <boost/math/constants/constants.hpp> //PI
 
//Ubitrack
#include "Numbers.h"
#include <utMath/Quaternion.h>

namespace Ubitrack { namespace Math { namespace Random {

/**
 * Function that generates a uniformly distributed quaternion.
 *
 * The functions implements the explanation regarding random unit quaternions from the following webside:
 * http://planning.cs.uiuc.edu/node198.html
 */

template< typename T = double > 
struct Quaternion
	: public std::unary_function< void, Math::Quaternion >
{

	public :
		Quaternion( )
			: std::unary_function< void, Math::Quaternion >( )
			{ };
	
		const Math::Quaternion operator()( void ) const
		{
			const T x = distribute_uniform< T >( 0, 1 );
			const T y = distribute_uniform< T >( 0, 1 );
			const T z = distribute_uniform< T >( 0, 1 );
			
			const T rootx = std::sqrt( x );
			const T rootxinv = std::sqrt( 1 - x );
			const T piz2 = 2 * boost::math::constants::pi< T >() * z;
			const T piy2 = 2 * boost::math::constants::pi< T >() * y;
			
			return Math::Quaternion( rootxinv * std::sin( piy2 ), rootxinv * std::cos( piy2 ), rootx * std::sin( piz2 ), rootx * std::cos( piz2 ) );
		}
};

}}} // namespace Ubitrack::Math::Random

#endif  // __H__RANDOM_ROTATIONS_H__

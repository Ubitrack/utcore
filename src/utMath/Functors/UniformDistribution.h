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

// TODO: implement module start/stop mechanics

/**
 * @ingroup math
 * @file
 * Implements a functor for generating uniformly distributed vectors.
 * Idealy this functor is combined with std::generate_n
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#ifndef __UBITRACK_UNIFORM_DISTRIUBUTION_H_INCLUDED__
#define __UBITRACK_UNIFORM_DISTRIUBUTION_H_INCLUDED__ 

// std
#include <iterator>
#include <functional>

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/constants/constants.hpp>

// Ubitrack
#include <utMath/Vector.h>
#include <utMath/Quaternion.h>
#include "../RandomNumbers.h"


namespace Ubitrack { namespace Math { namespace Functors {

/**
 * Functor that generates a vector of uniformly distributed random numbers.
 * Depending on the constructor call you have two options to define the limits.
 * Either you specify the limits once for all dimensions or you specify the 
 * limits for each dimension seperately.
 */
template< std::size_t N, typename T >
struct uniform_distribution
	: public std::unary_function< void, Math::Vector< N, T > >
{
	protected:
		const Math::Vector< N, T > m_min_range;
		const Math::Vector< N, T > m_max_range;
		
	public :
		uniform_distribution( const T min_range , const T max_range )
			: std::unary_function< void, Math::Vector< N, T > >( )
			, m_min_range( boost::numeric::ublas::scalar_vector< T >( N, std::min( min_range, max_range ) ) )
			, m_max_range( boost::numeric::ublas::scalar_vector< T >( N, std::max( min_range, max_range ) ) )
			{ };
			
		uniform_distribution( const Math::Vector< N, T > &min_range, const Math::Vector< N, T > &max_range )
			: std::unary_function< void, Math::Vector< N, T > >( )
			, m_min_range( min_range )
			, m_max_range( max_range )
			{ };

		const Math::Vector< N, T > operator()( void ) const
		{
			Math::Vector< N, T > vec; 
			for( std::size_t n( 0 ); n < N; ++n )
				vec( n ) = Math::distribute_uniform< T >( m_min_range( n ), m_max_range( n ) );
			return vec;
		}
};


/**
 * Function that generates a uniformly distributed quaternion.
 * the functions follows the explanation regarding random unit quaternions from the following webside:
 * http://planning.cs.uiuc.edu/node198.html
 */
template< typename T = double >
struct uniform_quaternion
	: public std::unary_function< void, Math::Quaternion >
{

	public :
		uniform_quaternion( )
			: std::unary_function< void, Math::Quaternion >( )
			{ };
	
		const Math::Quaternion operator()( void ) const
		{
			const T x = Math::distribute_uniform< T >( 0, 1 );
			const T y = Math::distribute_uniform< T >( 0, 1 );
			const T z = Math::distribute_uniform< T >( 0, 1 );
			
			const T rootx = std::sqrt( x );
			const T rootxinv = std::sqrt( 1 - x );
			const T piz2 = 2 * boost::math::constants::pi< T >() * z;
			const T piy2 = 2 * boost::math::constants::pi< T >() * y;
			
			return Math::Quaternion( rootxinv * std::sin( piy2 ), rootxinv * std::cos( piy2 ), rootx * std::sin( piz2 ), rootx * std::cos( piz2 ) );
		}
};

} } } //Ubitrack::Math::Functors

#endif // __UBITRACK_UNIFORM_DISTRIUBUTION_H_INCLUDED__
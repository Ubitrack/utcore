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
 * Implements functors for generating randomly distributed poses.
 *
 * The type of distribution can be either normal or uniform.
 * Ideally these functors are combined with \c std::generate_n.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#ifndef __UBITRACK_RANDOM_POSE_H_INCLUDED__
#define __UBITRACK_RANDOM_POSE_H_INCLUDED__ 

// std
#include <functional>


// Ubitrack
#include "Scalar.h"
#include "Rotation.h"
#include <utMath/Pose.h>
#include <utMath/Vector.h>


namespace Ubitrack { namespace Math { namespace Random {

/**
 * @ingroup math
 * Functor to generate randomly distributed poses.
 * @todo add functionality to produce normally distributed random poses.
 */
template< typename T >
struct Pose
{
	// /**
	 // * Functor that generates poses of normally distributed random numbers.
	 // *
	 // * Depending on the constructor call you have two options to define the limits.
	 // * Either you specify the mean value and standard deviation once for all dimensions
	 // * or your specify the mean value and standard deviation for each dimension separately.
	 // */
	// struct Normal
		// : public std::unary_function< void, Math::Pose >
	// {
		// protected:
			// const Math::Vector< N, T > m_mu;
			// const Math::Vector< N, T > m_sigma;
			
		// public :
			// Normal( const T mu, const T sigma )
				// : std::unary_function< void, Math::Vector< N, T > >( )
				// , m_mu( boost::numeric::ublas::scalar_vector< T >( N, mu ) )
				// , m_sigma( boost::numeric::ublas::scalar_vector< T >( N, sigma ) )
				// { };
		
			// Normal( const Math::Vector< N, T >& mu, const Math::Vector< N, T >& sigma )
				// : std::unary_function< void, Math::Vector< N, T > >( )
				// , m_mu( mu )
				// , m_sigma( sigma )
				// { };

			// const Math::Vector< N, T > operator()( void ) const
			// {
				// Math::Vector< N, T > vec; 
				// for( std::size_t n( 0 ); n < (N); ++n )
					// vec( n ) = distribute_normal< T >( m_mu( n ), m_sigma( n ) );
				// return Math::Vector< N, T >( vec );
			// }
	// };


	/**
	 * Functor that generates a uniformly distributed random pose.
	 *
	 * Depending on the constructor call you have two options to define the limits of the translation.
	 * Either you specify the limits once for all dimensions at the same time or you specify the 
	 * limits for each dimension separately.
	 */
	struct Uniform
		: public std::unary_function< void, Math::Pose >
	{
		protected:
			Random::Quaternion< double > randRotations;
			const Math::Vector< T, 3 > m_min_range;
			const Math::Vector< T, 3 > m_max_range;
			
		public :
			Uniform( const T min_range , const T max_range )
				: std::unary_function< void, Math::Pose >( )
				, m_min_range( boost::numeric::ublas::scalar_vector< T >( 3, std::min( min_range, max_range ) ) )
				, m_max_range( boost::numeric::ublas::scalar_vector< T >( 3, std::max( min_range, max_range ) ) )
				{ };
				
			Uniform( const Math::Vector< T, 3 > &min_range, const Math::Vector< T, 3 > &max_range )
				: std::unary_function< void, Math::Pose >( )
				, m_min_range( min_range )
				, m_max_range( max_range )
				{ };

			const Math::Pose operator()( void ) const
			{
				Math::Vector< T, 3 > vec; 
				vec( 0 ) = distribute_uniform< T >( m_min_range( 0 ), m_max_range( 0 ) );
				vec( 1 ) = distribute_uniform< T >( m_min_range( 1 ), m_max_range( 1 ) );
				vec( 2 ) = distribute_uniform< T >( m_min_range( 2 ), m_max_range( 2 ) );
				return Math::Pose( randRotations(), vec );
			}
	};
};

}}} //Ubitrack::Math::Random

#endif //__UBITRACK_RANDOM_POSE_H_INCLUDED__
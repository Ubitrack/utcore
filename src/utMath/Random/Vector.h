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
 * Implements functors for generating randomly distributed vectors.
 *
 * The type of distribution can be either normal or uniform.
 * Ideally this functor is combined with \c std::generate_n.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#ifndef __UBITRACK_RANDOM_VECTOR_H_INCLUDED__
#define __UBITRACK_RANDOM_VECTOR_H_INCLUDED__ 

// std
#include <functional>

// boost
#include <utMath/Vector.h>
#include <boost/numeric/bindings/lapack/posv.hpp>

// Ubitrack
#include "Scalar.h"
#include <utMath/Vector.h>
#include <utMath/Matrix.h>


namespace Ubitrack { namespace Math { namespace Random {


/**
* @ingroup math
* Functor to generate vectors of a random distribution
*/
template< typename T, std::size_t N >
struct Vector
{
	/**
	 * Functor that generates a vector of normally distributed random numbers.
	 *
	 * Depending on the constructor call you have two options to define the limits.
	 * Either you specify the mean value and standard deviation once for all dimensions
	 * or your specify the mean value and standard deviation for each dimension separately.
	 */
	struct Normal
		: public std::unary_function< void, Math::Vector< T, N > >
	{
		protected:
			const Math::Vector< T, N > m_mu;
			const Math::Vector< T, N > m_sigma;
			// 2013-11-07 CW: Commented out the drawing number from a covariance 
			// the intention was good, the realization was wrong
			// actually there is already a function to draw from a normal distribution
			// this function is used now.
			//maybe some time later I will take here once more...
			// Math::Matrix< T, N, N > m_sigma;

			// void prepareSigma()
			// {
				// //performs an inversion via cholesky decomposition
				// boost::numeric::bindings::lapack::potrf( 'L', m_sigma );

				// //set upper triangle of matrix to lower triangle
				// for( std::size_t r( 0 ); r < (N); ++r )
					// for( std::size_t c ( r + 1 ); c < (N); ++c )
						// m_sigma( r, c ) = m_sigma( c, r );
			// }
			
		public :
			Normal( const T mu, const T sigma )
				: std::unary_function< void, Math::Vector< T, N > >( )
				, m_mu( boost::numeric::ublas::scalar_vector< T >( N, mu ) )
				, m_sigma( boost::numeric::ublas::scalar_vector< T >( N, sigma ) )
				// , m_sigma( boost::numeric::ublas::identity_matrix< T >( N ) * sigma )
				{ 
					// prepareSigma();
				};
		
			Normal( const Math::Vector< T, N >& mu, const Math::Vector< T, N >& sigma )
				: std::unary_function< void, Math::Vector< T, N > >( )
				, m_mu( mu )
				, m_sigma( sigma )
				// , m_sigma( boost::numeric::ublas::zero_matrix< T >( N ) )
				{
					// for( std::size_t n( 0 ); n < N; ++n )
						// m_sigma( n, n ) = sigma( n );
					// prepareSigma();
				};

			Normal( const Math::Vector< T, N >& mu, const Math::Matrix< T, N, N >& sigma )
				: std::unary_function< void, Math::Vector< T, N > >( )
				, m_mu( mu )
				// , m_sigma( sigma )
				{
					///@TODO write correct algorithm how to draw from a random number from a given covariance
					//right now this initialization would be wrong...
					// prepareSigma();
				};

			const Math::Vector< T, N > operator()( void ) const
			{
				Math::Vector< T, N > vec; 
				for( std::size_t n( 0 ); n < N; ++n )
					vec( n ) = distribute_normal< T >( m_mu[ n ], m_sigma[ n ] );
				return vec;
				// return ( m_mu + boost::numeric::ublas::prod( m_sigma, vec ) );
				/*Math::Vector< T, N > vec; 
				for( std::size_t n( 0 ); n < (N); ++n )
					vec( n ) = distribute_normal< T >( m_mu( n ), m_sigma( n ) );
				return Math::Vector< T, N >( vec );*/
			}
	};


	/**
	 * Functor that generates a vector of uniformly distributed random numbers.
	 *
	 * Depending on the constructor call you have two options to define the limits.
	 * Either you specify the limits once for all dimensions or you specify the 
	 * limits for each dimension separately.
	 */

	struct Uniform
		: public std::unary_function< void, Math::Vector< T, N > >
	{
		protected:
			const Math::Vector< T, N > m_min_range;
			const Math::Vector< T, N > m_max_range;
			
		public :
			Uniform( const T min_range , const T max_range )
				: std::unary_function< void, Math::Vector< T, N > >( )
				, m_min_range( boost::numeric::ublas::scalar_vector< T >( N, std::min( min_range, max_range ) ) )
				, m_max_range( boost::numeric::ublas::scalar_vector< T >( N, std::max( min_range, max_range ) ) )
				{ };
				
			Uniform( const Math::Vector< T, N > &min_range, const Math::Vector< T, N > &max_range )
				: std::unary_function< void, Math::Vector< T, N > >( )
				, m_min_range( min_range )
				, m_max_range( max_range )
				{ };

			const Math::Vector< T, N > operator()( void ) const
			{
				Math::Vector< T, N > vec; 
				for( std::size_t n( 0 ); n < N; ++n )
					vec( n ) = distribute_uniform< T >( m_min_range( n ), m_max_range( n ) );
				return Math::Vector< T, N >( vec );
			}
	};
};

/** 
 * Function that produces a N-dimensional random vector of a given normal distribution.
 *
 * @tparam T type of distribution ( e.g. \c double or \float )
 * @tparam N dimension of the vector to be returned
 * @param mu mean value of normal distribution
 * @param sigma standard deviation of normal distribution
 * @return the random vector from the normal distribution
 */
template< typename T, std::size_t N >
Math::Vector< T, N > distribute_normal( const T mu , const T sigma )
{
	return typename Random::Vector< T, N >::Normal( mu, sigma )();
}

/** 
 * Function that produces a N-dimensional random vector of a given normal distribution.
 *
 * @tparam T type of distribution ( e.g. \c double or \float )
 * @tparam N dimension of the vector to be returned
 * @param mu mean value of normal distribution
 * @param sigma standard deviation of normal distribution
 * @return the random vector from the normal distribution
 */
template< typename T, std::size_t N >
Math::Vector< T, N > distribute_normal( const Math::Vector< T, N > &mu , const Math::Matrix< T, N, N >& sigma )
{
	return typename Random::Vector< T, N >::Normal( mu, sigma )();
}

/** 
 * Function that produces a N-dimensional random vector of a given uniform distribution.
 *
 * @tparam T type of distribution ( e.g. \c double or \float )
 * @tparam N dimension of the vector to be returned
 * @param min lower bound of uniform distribution
 * @param max upper bound of uniform distribution
 * @return the random vector with entries between min and max
 */
template< typename T, std::size_t N >
Math::Vector< T, N > distribute_uniform( const T min, const T max )
{
	return typename Random::Vector< T, N >::Uniform( min, max )();
}

}}} //Ubitrack::Math::Random

#endif //__UBITRACK_RANDOM_VECTOR_H_INCLUDED__

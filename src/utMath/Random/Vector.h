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
#include <boost/numeric/ublas/vector.hpp>
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
template< std::size_t N, typename T >
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
		: public std::unary_function< void, Math::Vector< N, T > >
	{
		protected:
			const Math::Vector< N, T > m_mu;
			Math::Matrix< N, N, T > m_sigma;

			void prepareSigma()
			{
				//performs an inversion via cholesky decomposition
				boost::numeric::bindings::lapack::potrf( 'L', m_sigma );

				//set upper triangle of matrix to zero
				for( std::size_t r( 0 ); r < (N); ++r )
					for( std::size_t c ( r + 1 ); c < (N); ++c )
						m_sigma( r, c ) = 0;
			}
			
		public :
			Normal( const T mu, const T sigma )
				: std::unary_function< void, Math::Vector< N, T > >( )
				, m_mu( boost::numeric::ublas::scalar_vector< T >( N, mu ) )
				, m_sigma( boost::numeric::ublas::identity_matrix< T >( N ) * sigma )
				{ 
					prepareSigma();
				};
		
			Normal( const Math::Vector< N, T >& mu, const Math::Vector< N, T >& sigma )
				: std::unary_function< void, Math::Vector< N, T > >( )
				, m_mu( mu )
				, m_sigma( boost::numeric::ublas::zero_matrix< T >( N ) )
				{
					for( std::size_t n( 0 ); n < N; ++n )
						m_sigma( n, n ) = sigma( n );
					prepareSigma();
				};

			Normal( const Math::Vector< N, T >& mu, const Math::Matrix< N, N, T >& sigma )
				: std::unary_function< void, Math::Vector< N, T > >( )
				, m_mu( mu )
				, m_sigma( sigma )
				{
					prepareSigma();
				};

			const Math::Vector< N, T > operator()( void ) const
			{
				Math::Vector< N, T > vec; 
				for( std::size_t n( 0 ); n < N; ++n )
					vec( n ) = distribute_normal< T >( 0, 1 );
				return ( m_mu + boost::numeric::ublas::prod( m_sigma, vec ) );
				/*Math::Vector< N, T > vec; 
				for( std::size_t n( 0 ); n < (N); ++n )
					vec( n ) = distribute_normal< T >( m_mu( n ), m_sigma( n ) );
				return Math::Vector< N, T >( vec );*/
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
		: public std::unary_function< void, Math::Vector< N, T > >
	{
		protected:
			const Math::Vector< N, T > m_min_range;
			const Math::Vector< N, T > m_max_range;
			
		public :
			Uniform( const T min_range , const T max_range )
				: std::unary_function< void, Math::Vector< N, T > >( )
				, m_min_range( boost::numeric::ublas::scalar_vector< T >( N, std::min( min_range, max_range ) ) )
				, m_max_range( boost::numeric::ublas::scalar_vector< T >( N, std::max( min_range, max_range ) ) )
				{ };
				
			Uniform( const Math::Vector< N, T > &min_range, const Math::Vector< N, T > &max_range )
				: std::unary_function< void, Math::Vector< N, T > >( )
				, m_min_range( min_range )
				, m_max_range( max_range )
				{ };

			const Math::Vector< N, T > operator()( void ) const
			{
				Math::Vector< N, T > vec; 
				for( std::size_t n( 0 ); n < N; ++n )
					vec( n ) = distribute_uniform< T >( m_min_range( n ), m_max_range( n ) );
				return Math::Vector< N, T >( vec );
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
Math::Vector< N, T > distribute_normal( const T mu , const T sigma )
{
	return typename Random::Vector< N, T >::Normal( mu, sigma )();
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
Math::Vector< N, T > distribute_normal( const Math::Vector< N, T > &mu , const Math::Matrix< N, N, T>& sigma )
{
	return typename Random::Vector< N, T >::Normal( mu, sigma )();
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
Math::Vector< N, T > distribute_uniform( const T min, const T max )
{
	return typename Random::Vector< N, T >::Uniform( min, max )();
}

}}} //Ubitrack::Math::Random

#endif //__UBITRACK_RANDOM_VECTOR_H_INCLUDED__

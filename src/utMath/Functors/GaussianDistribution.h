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
 * Implements a functor for generating normal distributed vectors.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#ifndef __UBITRACK_GAUSSIAN_DISTRIUBUTION_H_INCLUDED__
#define __UBITRACK_GAUSSIAN_DISTRIUBUTION_H_INCLUDED__ 

// std
#include <iterator>
#include <functional>

// boost
#include <boost/numeric/ublas/vector.hpp>

// Ubitrack
#include <utMath/Vector.h>
#include <utMath/Quaternion.h>
#include "../RandomNumbers.h"


namespace Ubitrack { namespace Math { namespace Functors {

/**
 * Functor that generates a vector of normally distributed random numbers.
 * Depending on the constructor call you have two options to define the limits.
 * Either you specify the mean value and standard deviation once for all dimensions
 * or your specify the mean value and standard deviation for each dimension seperately.
 */
template< std::size_t N, typename T >
struct gaussian_distribution
	: public std::unary_function< void, Math::Vector< N, T > >
{
	protected:
		const Math::Vector< N, T > m_mu;
		const Math::Vector< N, T > m_sigma;
		
	public :
		gaussian_distribution( const T mu, const T sigma )
			: std::unary_function< void, Math::Vector< N, T > >( )
			, m_mu( boost::numeric::ublas::scalar_vector< T >( N, mu ) )
			, m_sigma( boost::numeric::ublas::scalar_vector< T >( N, sigma ) )
			{ };
	
		gaussian_distribution( const Math::Vector< N, T > &mu, const Math::Vector< N, T > &sigma )
			: std::unary_function< void, Math::Vector< N, T > >( )
			, m_mu( mu )
			, m_sigma( sigma )
			{ };

		const Math::Vector< N, T > operator()( void ) const
		{
			Math::Vector< N, T > vec; 
			for( std::size_t n( 0 ); n < (N); ++n )
				vec( n ) = distribute_normal< T >( m_mu( n ), m_sigma( n ) );
			return Math::Vector< N, T >( vec );
		}
};

} } } //namespace Ubitrack::Math::Functors

#endif //__UBITRACK_GAUSSIAN_DISTRIUBUTION_H_INCLUDED__
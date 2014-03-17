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

#ifndef __UBITRACK_MATH_RANDOM_MATRIX_H_INCLUDED__
#define __UBITRACK_MATH_RANDOM_MATRIX_H_INCLUDED__ 

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
template< typename T, std::size_t M, std::size_t N >
struct Matrix
{

	/**
	 * Functor that generates a vector of uniformly distributed random numbers.
	 *
	 * Depending on the constructor call you have two options to define the limits.
	 * Either you specify the limits once for all dimensions or you specify the 
	 * limits for each dimension separately.
	 */

	struct Uniform
		: public std::unary_function< void, Math::Matrix< T, M, N > >
	{
		protected:
			const T m_min_range;
			const T m_max_range;
			
		public :
			Uniform( const T min_range , const T max_range )
				: std::unary_function< void, Math::Matrix< T, M, N > >( )
				, m_min_range( std::min( min_range, max_range ) ) 
				, m_max_range( std::max( min_range, max_range ) )
				{ };

			const Math::Matrix< T, M, N > operator()( void ) const
			{
				Math::Matrix< T, M, N > mat; 
				for( std::size_t m( 0 ); m < M; ++m )
					for( std::size_t n( 0 ); n < N; ++n )
						mat( m, n ) = distribute_uniform< T >( m_min_range, m_max_range );
				return Math::Matrix< T, M, N >( mat );
			}
	};
};

}}} //Ubitrack::Math::Random

#endif //__UBITRACK_MATH_RANDOM_MATRIX_H_INCLUDED__

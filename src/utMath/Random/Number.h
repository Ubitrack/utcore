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
 * generates a one-dimensional random number
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

 
 
#ifndef __H__RANDOM_NUMBERS_H__
#define __H__RANDOM_NUMBERS_H__

// boost random numbers available from version 1.31.xx on
#if (BOOST_VERSION / 100000) > 0 && (BOOST_VERSION / 100 % 1000 ) > 30
	#define RANDOM_BOOST
	#include <boost/random/uniform_01.hpp>
	#include <boost/random/uniform_real.hpp>
	#include <boost/random/variate_generator.hpp>
	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/normal_distribution.hpp>
#else
///@todo C++0x also offers random variables, add this as a possibility as well
#include <cstdlib>
#endif

namespace Ubitrack { namespace Math { namespace Random {

#ifdef RANDOM_BOOST
static boost::mt19937 RNG; // Random Number Generator
#endif

/** 
 * Function that produces a one dimensional random number of a given normal distribution
 *
 * @param mu mean value of normal distribution
 * @param sigma standard deviation of normal distribution
 * @return the random number
 */
template< typename T >
T distribute_normal( const T mu , const T sigma )
{
#ifdef RANDOM_BOOST
	boost::normal_distribution< T > normalDist( mu, sigma ); //Normal (=Gaussian) distribution
	boost::variate_generator< boost::mt19937&, boost::normal_distribution< T > > Generator( RNG, normalDist );
	return Generator();
#else
	/// @todo: change this, right now it's no good idea for normal distribution
	return rand()
#endif
}

/** 
 *
 * Function that produces a one dimensional random number of a given uniform distribution
 *
 * @tparam T type of distribution ( e.g. \c double or \float )
 * @param min lower bound of uniform distribution
 * @param max upper bound of uniform distribution
 * @return the random number between min and max
 */
template< typename T >
T distribute_uniform( const T min, const T max )
{
#ifdef RANDOM_BOOST
	boost::uniform_real< T > uniformDist( min, max ); // Uniform distribution
	boost::variate_generator < boost::mt19937&, boost::uniform_real< T > > Generator( RNG, uniformDist );
	return Generator();
#else
	/// @todo: change this, right now it's no good idea for uniform distribution
	return rand()% max + min;
#endif
}

}}} // namespace Ubitrack::Math::Random

#endif  // __H__RANDOM_NUMBERS_H__

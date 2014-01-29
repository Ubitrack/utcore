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


#ifndef __UBITRACK_MATH_STOCHASTIC_GAUSSIAN_H__
#define __UBITRACK_MATH_STOCHASTIC_GAUSSIAN_H__
 
// Ubitrack
#include <utUtil/Exception.h>
// #include <utUtil/StaticAssert.h>

//std 
#include <iomanip> // stream output

namespace Ubitrack{ namespace Math { namespace Stochastic {

/// @internal a generic Gaussian distribution type
template< typename T, std::size_t N >
struct Gaussian
{
	/** simple typedef to type of Gaussian structure. */
	typedef T value_type;
	typedef std::size_t size_type;
	static const size_type size = N; 
		
	/** the mean value */
	T mean[ N ];
	
	/** the covariance matrix, upper triangle equals the lower triangle */
	T covariance [ N * N ];
	
	/** the root of the sum of squared diagonal entries of the covariance */
	T variance;
	
	/** the sum of squared diagonal entries of the covariance */
	T squarredVariance;
};


/// @internal a function that estimates a Gaussian distribution
template< typename T, std::size_t N, typename InputIterator1, typename InputIterator2 >
bool estimate_gaussian( const InputIterator1 itBegin, const InputIterator1 itEnd, const InputIterator2 itWeights, Gaussian< T, N > &gaussian )
{
	const std::size_t n = std::distance( itBegin, itEnd );
	
	if( n == 0 )// no value exit early
		return false;
		
	// set the covariance to zero elements
	gaussian.variance = gaussian.squarredVariance = 0;
	std::fill( gaussian.covariance, gaussian.covariance+(N*N), static_cast< T >( 0 ) );
	if( n == 1 ) //only one value, simple case
	{
		for( std::size_t i( 0 ); i < N; ++i )
			gaussian.mean[ i ] = (*itBegin)[ i ];
		return true;
	}
	// since we got that far also zero out the mean value
	std::fill( gaussian.mean, gaussian.mean+N, static_cast< T >( 0 ) );
	// Calculate (weighted) mean value
	{
		InputIterator2 weightIter( itWeights );
		for( InputIterator1 valIter ( itBegin ); valIter != itEnd; ++valIter, ++weightIter )
			for( std::size_t i( 0 ); i < N; ++i )
				gaussian.mean[ i ] += (*weightIter) * (*valIter)[ i ];
	}
		
	// Calculate (weighted) covariance
	{
		InputIterator2 weightIter( itWeights );
		for( InputIterator1 valIter ( itBegin ); valIter != itEnd; ++valIter, ++weightIter )
			for( std::size_t i1( 0 ); i1 < N; ++i1 )
				for( std::size_t i2( i1 ); i2 < N; ++i2 )
					gaussian.covariance[ i2*N+i1 ] = gaussian.covariance[ i1*N+i2 ] += ( *weightIter ) * ((*valIter)[ i1 ] - gaussian.mean[ i1 ] ) * ( (*valIter)[ i2 ] - gaussian.mean[ i2 ] );
	}
	
	{
		//sum up the squared diagonal entries
		for( std::size_t i( 0 ); i < N; ++i )	
			gaussian.squarredVariance += ( gaussian.covariance[i*N+i] * gaussian.covariance[i*N+i] );
			
		gaussian.variance = std::sqrt( gaussian.squarredVariance );
	}
	return true;
};

template< typename T, std::size_t N, typename InputIterator >
bool estimate_gaussian( const InputIterator itBegin, const InputIterator itEnd, Gaussian< T, N > &gaussian )
{
	const std::size_t n = std::distance( itBegin, itEnd );
		
	if( n == 0 )// no value exit early
		return;

	// set the covariance to zero elements
	gaussian.variance = gaussian.squarredVariance = 0;
	std::fill( gaussian.covariance, gaussian.covariance+(N*N), static_cast< T >( 0 ) );
	if( n == 1 ) //only one value, simple case
	{
		for( std::size_t i( 0 ); i < N; ++i )
			gaussian.mean[ i ] = (*itBegin)[ i ];
		return true;
	}
	// since we got that far also zero out the mean value
	std::fill( gaussian.mean, gaussian.mean+N, static_cast< T >( 0 ) );
	
	// calculate mean value
	{
		for( InputIterator it ( itBegin ); it != itEnd; ++it )
			for( std::size_t i( 0 ); i < N; ++i )
				gaussian.mean[ i ] += (*it)[ i ];
	
		for( std::size_t i( 0 ); i < N; ++i )
			gaussian.mean[ i ] /= n;
	}
	
	// calculate covariance
	{
		for( InputIterator it ( itBegin ); it != itEnd; ++it )
			for( std::size_t i1( 0 ); i1 < N; ++i1 )
				for( std::size_t i2( i1 ); i2 < N; ++i2 )
					gaussian.covariance[ i2*N+i1 ] = gaussian.covariance[ i1*N+i2 ] += ((*it)[ i1 ] - gaussian.mean[ i1 ]) * ( (*it)[ i2 ] - gaussian.mean[ i2 ]);
				
		for( std::size_t i( 0 ); i < (N*N); ++i )
			gaussian.covariance[ i ] /= n;
	}
	
	{
		//sum up the squared diagonal entries
		for( std::size_t i( 0 ); i < N; ++i )	
			gaussian.squarredVariance += ( gaussian.covariance[i*N+i] * gaussian.covariance[i*N+i] );
			
		gaussian.variance = std::sqrt( gaussian.squarredVariance );
	}
	return true;
	
	
	// calculate mean value + covariance
	// {
		
		// for( FwdValueIter valIter ( itValBegin ); valIter != itValEnd; ++valIter )
			// for( std::size_t i1( 0 ); i1 < N; ++i1 )
			// {
				// mean[ i1 ] += (*weightIter) * (*valIter)[ i1 ];
				// for( std::size_t i2( i1 ); i2 < N; ++i2 )
					// covariance[ i2*N+i1 ] = covariance[ i1*N+i2 ] += * ((*valIter)[ i1 ]  * ( (*valIter)[ i2 ] ));
			// }
		
		// normalizeMean( n );
		// normalizeCovariance( n );
		
		// for( std::size_t i1( 0 ); i1 < N; ++i1 )
			// for( std::size_t i2( 0 ); i2 < N; ++i2 )
				// covariance[ i1*N+i2 ] -= mean[ i1 ] * mean[ i2 ];
	// }
};


template< typename T, std::size_t N, typename InputIterator1, typename InputIterator2 >
bool estimate_gaussian( const InputIterator1 itBegin, const InputIterator1 itEnd, const InputIterator2 itIndices, const typename std::iterator_traits< InputIterator2 >::value_type comp_value, Gaussian< T, N > &gaussian )
{
	const std::size_t n = std::distance( itBegin, itEnd );
		
	if( n == 0 )// no value exit early
		return false;

	// set the covariance to zero elements
	gaussian.variance = gaussian.squarredVariance = 0;
	std::fill( gaussian.covariance, gaussian.covariance+(N*N), static_cast< T >( 0 ) );
	if( n == 1 ) //only one value, simple case
	{
		for( std::size_t i( 0 ); i < N; ++i )
			gaussian.mean[ i ] = (*itBegin)[ i ];
		return true;
	}
	// since we got that far also zero out the mean value
	std::fill( gaussian.mean, gaussian.mean+N, static_cast< T >( 0 ) );
	
	std::size_t count( 0 );
	// calculate mean value
	{
		InputIterator2 indexIter( itIndices );
		for( InputIterator1 it ( itBegin ); it != itEnd; ++it, ++indexIter )
			if( (*indexIter) == comp_value )
			{
				count++;
				for( std::size_t i( 0 ); i < N; ++i )
					gaussian.mean[ i ] += (*it)[ i ];
			}
		for( std::size_t i( 0 ); i < N; ++i )
			gaussian.mean[ i ] /= count;
	}
	
	// calculate covariance
	{
		InputIterator2 indexIter( itIndices );
		for( InputIterator1 it ( itBegin ); it != itEnd; ++it, ++indexIter )
			if( (*indexIter) == comp_value )
				for( std::size_t i1( 0 ); i1 < N; ++i1 )
					for( std::size_t i2( i1 ); i2 < N; ++i2 )
						gaussian.covariance[ i2*N+i1 ] = gaussian.covariance[ i1*N+i2 ] += ((*it)[ i1 ] - gaussian.mean[ i1 ]) * ( (*it)[ i2 ] - gaussian.mean[ i2 ]);
				
		for( std::size_t i( 0 ); i < (N*N); ++i )
			gaussian.covariance[ i ] /= count;
	}
	
	{
		//sum up the squared diagonal entries
		for( std::size_t i( 0 ); i < N; ++i )	
			gaussian.squarredVariance += ( gaussian.covariance[i*N+i] * gaussian.covariance[i*N+i] );
			
		gaussian.variance = std::sqrt( gaussian.squarredVariance );
	}
	return true;
}



/** @internal overrides the stream output to have nicely aligned data */
template< typename T, std::size_t N >
std::ostream& operator<<( std::ostream& s, const Gaussian< T, N >& gauss )
{
	// s << gauss.variance << "\n" << gauss.squarredVariance << "\n";
	for( std::size_t i1( 0 ); i1<N; ++i1 )
	{
		s << std::setfill(' ')
		<< std::setw(10)
		<< std::fixed
		<< std::setprecision(2)
		<< gauss.mean[ i1 ] << " [ ";
		for( std::size_t i2( 0 ); i2<N; ++i2 )
			s << std::setw( 10 )
			<< std::fixed
			<< std::setprecision(4)
			<< gauss.covariance[ i1*N+i2 ] ;
		s << " ]" << std::endl;
	}
	return s;
}


} } } // namespace Ubitrack::Math::Stochastic

#endif //__UBITRACK_MATH_STOCHASTIC_GAUSSIAN_H__
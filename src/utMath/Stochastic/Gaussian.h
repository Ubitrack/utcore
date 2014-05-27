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
	/** typedef to built-in type of Gaussian distribution */
	typedef T value_type;
	
	/** typedef to dimension type of Gaussian distribution */
	typedef std::size_t size_type;
	
	/** dimension of the Gaussian distribution */
	static const size_type size = N; 
		
	/** the mean value */
	T mean[ N ];
	
	/** the covariance matrix, upper triangle equals the lower triangle */
	T covariance [ N * N ];
	
	/** the root of the sum of squared diagonal entries of the covariance */
	T variance;
	
	/** the sum of squared diagonal entries of the covariance */
	T squaredVariance;
};

/// @internal a function that sets all values of the gaussian to zero
template< typename T, std::size_t N >
void reset( Gaussian< T, N > &gaussian )
{
	// set all values of the gaussian type to zero
	std::fill( gaussian.mean, gaussian.mean+N, static_cast< T >( 0 ) );
	std::fill( gaussian.covariance, gaussian.covariance+(N*N), static_cast< T >( 0 ) );
	gaussian.variance = gaussian.squaredVariance = 0;
}


/// @internal a function that estimates a (weighted) Gaussian distribution
template< typename T, std::size_t N, typename InputIterator1, typename InputIterator2 >
bool estimate_gaussian( const InputIterator1 itBegin, const InputIterator1 itEnd, const InputIterator2 itWeights, Gaussian< T, N > &gaussian )
{
	// zero out the gaussian
	reset( gaussian );
	
	std::size_t n( 0 );
	{ // calculate (weighted) mean value
		InputIterator2 weightIter( itWeights );
		for( InputIterator1 valIter ( itBegin ); valIter != itEnd; ++valIter, ++weightIter, ++n )
			for( std::size_t i( 0 ); i < N; ++i )
				gaussian.mean[ i ] += (*weightIter) * (*valIter)[ i ];
				
		if( n < 1 )
			return false;
	}
		
	
	{ // calculate (weighted) covariance
		InputIterator2 weightIter( itWeights );
		for( InputIterator1 valIter ( itBegin ); valIter != itEnd; ++valIter, ++weightIter )
			for( std::size_t i1( 0 ); i1 < N; ++i1 )
				for( std::size_t i2( i1 ); i2 < N; ++i2 )
					gaussian.covariance[ i2*N+i1 ] = gaussian.covariance[ i1*N+i2 ] += ( *weightIter ) * ((*valIter)[ i1 ] - gaussian.mean[ i1 ] ) * ( (*valIter)[ i2 ] - gaussian.mean[ i2 ] );
	}
	
	{//sum up the squared diagonal entries
		for( std::size_t i( 0 ); i < N; ++i )	
			gaussian.squaredVariance += ( gaussian.covariance[i*N+i]  );
			
		gaussian.variance = std::sqrt( gaussian.squaredVariance );
	}
	return true;
};

/// @internal a function that estimates a Gaussian distribution
template< typename T, std::size_t N, typename InputIterator >
bool estimate_gaussian( const InputIterator itBegin, const InputIterator itEnd, Gaussian< T, N > &gaussian )
{
	reset( gaussian );
	
	std::size_t n( 0 );
	
	{ // calculate mean value
		for( InputIterator it ( itBegin ); it != itEnd; ++it, ++n )
		{
			std::cout << n << "\n";
			for( std::size_t i( 0 ); i < N; ++i )
				gaussian.mean[ i ] += (*it)[ i ];
		}
		if( n < 1 )
			return false;
			
		for( std::size_t i( 0 ); i < N; ++i )
			gaussian.mean[ i ] /= n;
	}
	
	
	{ // calculate covariance
		for( InputIterator it ( itBegin ); it != itEnd; ++it )
			for( std::size_t i1( 0 ); i1 < N; ++i1 )
				for( std::size_t i2( i1 ); i2 < N; ++i2 )
					gaussian.covariance[ i2*N+i1 ] = gaussian.covariance[ i1*N+i2 ] += ((*it)[ i1 ] - gaussian.mean[ i1 ]) * ( (*it)[ i2 ] - gaussian.mean[ i2 ]);
				
		for( std::size_t i( 0 ); i < (N*N); ++i )
			gaussian.covariance[ i ] /= n;
	}
	
	std::cout <<"Finsished With Variance" << std::endl;
	
	{ //sum up the squared diagonal entries
		for( std::size_t i( 0 ); i < N; ++i )	
			gaussian.squaredVariance += ( gaussian.covariance[i*N+i] );
			
		gaussian.variance = std::sqrt( gaussian.squaredVariance );
	}
	return true;
}

/// @internal a function that estimates a Gaussian distribution
template< typename T, std::size_t N, typename InputIterator >
bool estimate_gaussian_fast( const InputIterator itBegin, const InputIterator itEnd, Gaussian< T, N > &gaussian )
{
	reset( gaussian );
	std::size_t n( 0 );
	//calculate mean value + covariance
	// {
		
		// for( InputIterator valIter ( itBegin ); valIter != itEnd; ++valIter, ++n )
			// for( std::size_t i1( 0 ); i1 < N; ++i1 )
			// {
				// gaussian.mean[ i1 ] += (*valIter)[ i1 ];
				// for( std::size_t i2( i1 ); i2 < N; ++i2 )
					// gaussian.covariance[ i2*N+i1 ] = covariance[ i1*N+i2 ] += * ((*valIter)[ i1 ]  * ( (*valIter)[ i2 ] ));
			// }
		
		// normalizeMean( n );
		// normalizeCovariance( n );
		
		// for( std::size_t i1( 0 ); i1 < N; ++i1 )
			// for( std::size_t i2( 0 ); i2 < N; ++i2 )
				// covariance[ i1*N+i2 ] -= mean[ i1 ] * mean[ i2 ];
	// }
};


template< typename T, std::size_t N, typename InputIterator1, typename InputIterator2 >
bool estimate_gaussian_index( const InputIterator1 itBegin, const InputIterator1 itEnd, const InputIterator2 itIndices, const typename std::iterator_traits< InputIterator2 >::value_type comp_value, Gaussian< T, N > &gaussian )
{
	reset( gaussian );
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
			gaussian.squaredVariance += ( gaussian.covariance[i*N+i] );
			
		gaussian.variance = std::sqrt( gaussian.squaredVariance );
	}
	return true;
}



/** @internal overrides the stream output to have nicely aligned data */
template< typename T, std::size_t N >
std::ostream& operator<<( std::ostream& s, const Gaussian< T, N >& gauss )
{
	s << "Variance   : " << gauss.variance << "\nVariance^2 : " << gauss.squaredVariance << "\n";
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
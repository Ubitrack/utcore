/**
 * @ingroup math
 * @file
 * Functions for Correlation
 *
 * @author Manuel Huber <huberma@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#ifndef __UBITRACK_MATH_STOCHASTIC_CORRELATION_H_INCLUDED__
#define __UBITRACK_MATH_STOCHASTIC_CORRELATION_H_INCLUDED__

// std
#include <assert.h> //assert
#include <iterator> //std::iterator_traits
#include <cmath> // std::sqrt

namespace Ubitrack { namespace Math { namespace Stochastic {

/**
 * Correlation function that estimates the correlation factor of two one-dimensional container sequences.
 *
 * @tparam InputIterator type of iteterator pointing to both sequence containers including one-dimensional data for correlation
 * @param itBegin1 iterator pointing to first element in the first sequence container data.
 * @param itEnd1 iterator pointing to end of first sequence container data.
 * @param itBegin2 iterator pointing to first element in the second sequence container data.
 * @param itEnd2 iterator pointing to end of second sequence container data.
 * @return a value between [-1;1] that signs the correlation between the two data sequences.
 */
template< typename InputIterator >
typename std::iterator_traits< InputIterator >::value_type correlation( const InputIterator itBegin1, const InputIterator itEnd1
								, const InputIterator itBegin2 , const InputIterator itEnd2 )
{
	typedef typename std::iterator_traits< InputIterator >::value_type precision_type;
	
	const std::size_t n1 = std::distance( itBegin1, itEnd1 );
	const std::size_t n2 = std::distance( itBegin2, itEnd2 );
	
	if( (!n1) && (!n2) )
		return 1;
		
	assert( n1 ); //<- should be greater than zero
	assert( n2 ); //<- should be greater than zero
	
	// determine maximal length of chunks
	const std::size_t n = std::min( n1, n2 );
	
	// set the end of the first iterator to determined length
	InputIterator iEnd ( itBegin1 );
	std::advance( iEnd, n );
	const InputIterator cEnd ( iEnd );
	// const cEnd = (n == n1) ? itEnd1 : itBegin1+n; // <- no good idea, works only for random access iterators
	
	// determine mean value for the two chunks
	precision_type m1 = 0;
	precision_type m2 = 0;
	for( InputIterator iter1 = itBegin1, iter2 = itBegin2;  iter1 != cEnd ; ++iter1, ++iter2 )
	{
		m1 += (*iter1);
		m2 += (*iter2);
	}
	m1 /= n;
	m2 /= n;

	// determine the correlation factors
	precision_type res = 0;
	precision_type var1 = 0;
	precision_type var2 = 0;
	for( InputIterator iter1 = itBegin1, iter2 = itBegin2;  iter1 != cEnd ; ++iter1, ++iter2 )
	{
		const precision_type diff1 = (*iter1) - m1;
		const precision_type diff2 = (*iter2) - m2;
		
		res += diff1 * diff2;
		var1 += diff1 * diff1;
		var2 += diff2 * diff2;
	}

	const precision_type denom = std::sqrt( var1*var2 );

	return res/denom;
}

}}}; // namespace Ubitrack::Math::Stochastic


#endif // __UBITRACK_MATH_STOCHASTIC_CORRELATION_H_INCLUDED__

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
 * Ransac algorithm
 * @author Daniel Pustka <pustka@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */


#ifndef __UBITRACK_MATH_OPTIMIZATION_RANSAC_INCLUDED__
#define __UBITRACK_MATH_OPTIMIZATION_RANSAC_INCLUDED__

#include <utCore.h>
#include "Optimization.h"

#include <vector>
#include <stdlib.h>
#include <iterator> // std::iterator_traits
#include <algorithm> // std::generate_n
#include <stdexcept>


namespace Ubitrack { namespace Math { namespace Optimization {

/**
 * Parameter structure for the RANSAC algorithm framework
 *
 * This struct holds all parameters to set the appropriate parameters
 * of the ransac algorithm framework.
 */
template< typename T >
struct RansacParameter
{
public:
	/** type of all the unsigned integer values within this class */
	typedef std::size_t size_type;
	
	/** the floating point type used for numerical error expressions */
	typedef T value_type;

	const value_type threshold;
	const size_type setSize;
	const size_type nMinInlier;
	const size_type nMaxIterations;
	
	RansacParameter( const T fThreshold, const std::size_t n, const std::size_t minInlier, const std::size_t maxRuns )
		: threshold ( fThreshold )
		, setSize ( n )
		, nMinInlier ( minInlier )
		, nMaxIterations ( maxRuns )
		{};
	
	/**
	 * Constructor accepting intuitive parameters
	 *
	 *  @param fThreshold threhold to decide wether an value is accepeted tu support the hypothesis or not
	 *  @param nMinSet amount of values needed to estimate a solution by the given problem
	 *  @param n amount of values provided to estimate a solution for the given problem
	 *  @param percentOutlier signs the percentage of expected outlier ( ranges from 0 to 1 )
	 *  @param percentSucess probability that signs how intensive should be found a solution (affects number of iterations)
	 */
	RansacParameter( const T fThreshold, const std::size_t nMinSet, const std::size_t n, const T percentOutlier, const T percentSucess = 0.99 )
		: threshold ( fThreshold )
		, setSize ( nMinSet )
		, nMinInlier ( (1.-percentOutlier) * n  )
		, nMaxIterations ( static_cast< size_type >( 1+std::log( 1-percentSucess) / ( std::log( 1-std::pow( 1-percentOutlier, static_cast< int >( nMinSet ) ) ) ) ) )
		{};
};


/// @internal index generator helper struct
struct IndexGenerator
{
protected:
	std::size_t index;
	
public:
	IndexGenerator()
		: index( 0 )
	{}
	
	std::size_t operator()()
	{
		return index++;
	}
};

/**
 * RANSAC algorithm (for one-parameter problems)
 *
 * @tparam InputIterator describes the type of container iterator that points to the values
 * @tparam ResultType the result type of the solution formulation
 * @tparam T describes the numeric type used for error calculation (usually \c float or \c double )
 * @tparam RansacFunctor the type of the struct/class that should include \c Estimator and \c Evaluator functor object to estimate the solution and validate it
 * @param iBegin an \c iterator point to the first element of a container including the values
 * @param iEnd an \c iterator point to the final element of a container including the values
 * @param result returns the best estimated result for the given problem and parameter set
 * @param model an instance of the struct/class that includes the Estimator and Evaluator FunctorObjects that describe the solution of a problem and it's validation
 * @param params an instance of the object containing the algorithms parametrization
 * @return 0 (failure) or number of inlier on success
*/
template< class InputIterator, class ResultType, typename T, class RansacFunctor >
std::size_t ransac( const InputIterator iBegin, const InputIterator iEnd
	, ResultType& result
	, const RansacFunctor& model
	, const RansacParameter< T >& params )
{
	typedef typename std::iterator_traits< InputIterator >::value_type value_type;
	typedef typename value_type::value_type numeric_type;
	
	typedef std::vector< value_type > list_type;
	
	// estimate number of parameter list
	const std::size_t nValues = std::distance( iBegin, iEnd );
	assert( params.nMinInlier <= nValues );
	
	OPT_LOG_DEBUG( "RANSAC with " << nValues << " values , " << params.nMinInlier << " inlier required" );
	
	// includes every indices once, needed to shuffle for random order
	std::vector< std::size_t > iShuffle;
	iShuffle.reserve( nValues );
	std::generate_n( std::back_inserter( iShuffle ), nValues, IndexGenerator() );
	
	// set of pointer to inlier
	std::vector< InputIterator > iInliers;
	iInliers.reserve( nValues );

	// amount of best inlier so far
	std::size_t nBestInliers = 0;
	
	std::vector< InputIterator > iBestInliers;
	iBestInliers.reserve( nValues );
	
	std::size_t iRun ;
	for( iRun = 0; iRun < params.nMaxIterations; iRun++ )
	{
		OPT_LOG_TRACE( "RANSAC iteration " << iRun + 1 );
		
		
		// generate random set
		list_type list;
		list.reserve( params.setSize );
		
		// shuffle the indices little bit
		std::random_shuffle( iShuffle.begin(), iShuffle.end() );
		std::vector< std::size_t >::iterator itSelected = iShuffle.begin();
		std::vector< std::size_t >::iterator itEnd = iShuffle.begin();
		std::advance( itEnd, params.setSize );
		
		for ( ; itSelected < itEnd; ++itSelected )
		{
			InputIterator it ( iBegin );
			std::advance( it, (*itSelected) );
			list.push_back( *it );
		}
		
		// compute hypothesis
		ResultType hypothesis;
		if( ! typename RansacFunctor::Estimator()( hypothesis, list.begin(), list.end() ) )
		{
			OPT_LOG_TRACE( "fast forward, no estimation possible" );
			continue;
		}
		// count inlier
		std::size_t nInlier = 0;
		T fInlierDist = 0;
		iInliers.clear();
		
		InputIterator it ( iBegin );
		for ( std::size_t i = 0; i < nValues && int( nValues - i ) >= int( params.nMinInlier - nInlier ); i++, ++it )
		{
			const T d = typename RansacFunctor::Evaluator()( hypothesis, *it );
			if( d < params.threshold )
			{
				iInliers.push_back( it );
				nInlier++;
				fInlierDist += d;
			}
		}
		if( !nInlier )
		{
			OPT_LOG_TRACE( "fast forward, no inlier found" );
			continue;
		}	
		OPT_LOG_TRACE( nInlier << " inlier, avg dist=" << fInlierDist / nInlier );

		// save inlier set if we reached the required number or are better than a previous run
		if ( nInlier >= params.nMinInlier && nInlier > nBestInliers )
		{
			nBestInliers = nInlier;
			iBestInliers.swap( iInliers );
		}

		// stop after nMinRun iterations if the required number of inlier was found
		if ( nBestInliers >= params.nMinInlier )
			break;
	}

	if ( nBestInliers >= params.nMinInlier )
	{
		// compute final result
		list_type list;
		list.reserve( nBestInliers );
		
		typename std::vector< InputIterator >::const_iterator it = iBestInliers.begin();
		typename std::vector< InputIterator >::const_iterator itEnd = iBestInliers.end();
		for ( ; it < itEnd; ++it )
			list.push_back( *(*it) );

		typename RansacFunctor::Estimator()( result, list.begin(), list.end() );
		OPT_LOG_DEBUG( "Estimated " << nBestInliers << " inlier after " << iRun + 1 << " iterations."  );
		return nBestInliers;
	}
	
	OPT_LOG_DEBUG( "RANSAC: Not enough inlier found" );
	return 0;
}


/**
 * RANSAC algorithm (for two-parameter problems)
 *
 * @tparam InputIterator1 describes the type of container iterator that points to the first type of values
 * @tparam InputIterator2 describes the type of container iterator that points to the second type of values
 * @tparam ResultType the result type of the solution formulation
 * @tparam T describes the numeric type used for error calculation (usually \c float or \c double )
 * @tparam RansacFunctor the type of the struct/class that should include \c Estimator and \c Evaluator functor object to estimate the solution and validate it
 * @param iBegin1 an \c iterator that points to the first element of a container including the first problem values
 * @param iEnd1 an \c iterator that points to the final element of a container including the first problem values
 * @param iBegin2 an \c iterator that points to the first element of a container including the second problem values
 * @param iEnd2 an \c iterator that points to the final element of a container including the second problem values
 * @param result returns the best estimated result for the given problem and parameter set
 * @param model an instance of the struct/class that includes the Estimator and Evaluator FunctorObjects that describe the solution of a problem and it's validation
 * @param params an instance of the object containing the algorithms parametrization
 * @return 0 (failure) or number of inlier on success
*/
template< typename InputIterator1, typename InputIterator2, class ResultType, typename T, class RansacFunctor >
std::size_t ransac( const InputIterator1 iBegin1, const InputIterator1 iEnd1
		, const InputIterator2 iBegin2, const InputIterator2 iEnd2
		, ResultType& result
		, const RansacFunctor& model
		, const RansacParameter< T >& params )
{
	typedef typename std::iterator_traits< InputIterator1 >::value_type value_type;
	typedef typename value_type::value_type numeric_type;
	
	typedef std::vector< value_type > list_type;
	
	// estimate number of parameter list
	const std::size_t nValues = std::distance( iBegin1, iEnd1 );
	assert( params.nMinInlier <= nValues );
	
	OPT_LOG_DEBUG( "RANSAC with " << nValues << " values , " << params.nMinInlier << " inlier required" );
	
	// includes every indices once, needed to shuffle for random order
	std::vector< std::size_t > iShuffle;
	iShuffle.reserve( nValues );
	std::generate_n( std::back_inserter( iShuffle ), nValues, IndexGenerator() );
	
	// indices to inlier
	std::vector< std::size_t > iInliers;
	iInliers.reserve( nValues );

	// amount of best inlier so far
	std::size_t nBestInliers = 0;
	
	// indices to best set of inlier
	std::vector< std::size_t > iBestInliers;
	iBestInliers.reserve( nValues );
	
	std::size_t iRun ;
	for( iRun = 0; iRun < params.nMaxIterations; iRun++ )
	{
		OPT_LOG_TRACE( "RANSAC iteration " << iRun + 1 );
		
		
		// generate random set
		list_type list1;
		list_type list2;
		list1.reserve( params.setSize );
		list2.reserve( params.setSize );
		
		std::random_shuffle( iShuffle.begin(), iShuffle.end() );
		std::vector< std::size_t >::iterator itSelected = iShuffle.begin();
		std::vector< std::size_t >::iterator itEnd = iShuffle.begin();
		std::advance( itEnd, params.setSize );
		
		for ( ; itSelected < itEnd; ++itSelected )
		{
			InputIterator1 it1 ( iBegin1 );
			InputIterator2 it2 ( iBegin2 );
			std::advance( it1, (*itSelected) );
			std::advance( it2, (*itSelected) );
			list1.push_back( *it1 );
			list2.push_back( *it2 );
		}
		
		// compute hypothesis
		ResultType hypothesis;
		if( ! typename RansacFunctor::Estimator()( hypothesis, list1.begin(), list1.end(), list2.begin(), list2.end() ) )
		{
			OPT_LOG_TRACE( "fast forward, no estimation possible" );
			continue;
		}
		// count inlier
		std::size_t nInlier = 0;
		T fInlierDist = 0;
		iInliers.clear();
		
		InputIterator1 it1 ( iBegin1 );
		InputIterator2 it2 ( iBegin2 );
		for ( std::size_t i = 0; i < nValues && int( nValues - i ) >= int( params.nMinInlier - nInlier ); i++, ++it1, ++it2 )
		{
			const T d = typename RansacFunctor::Evaluator()( hypothesis, *it1, *it2 );
			if( d < params.threshold )
			{
				const std::size_t index = std::distance( iBegin1, it1 );
				iInliers.push_back( index );
				nInlier++;
				fInlierDist += d;
			}
		}
		if( !nInlier )
		{
			OPT_LOG_TRACE( "fast forward, no inlier found" );
			continue;
		}	
		OPT_LOG_TRACE( nInlier << " inlier, avg dist=" << fInlierDist / nInlier );

		// save inlier set if we reached the required number or are better than a previous run
		if ( nInlier >= params.nMinInlier && nInlier > nBestInliers )
		{
			nBestInliers = nInlier;
			iBestInliers.swap( iInliers );
		}

		// stop after nMinRun iterations if the required number of inlier was found
		if ( nBestInliers >= params.nMinInlier )
			break;
	}

	if ( nBestInliers >= params.nMinInlier )
	{
		// compute final result
		list_type list1;
		list_type list2;
		list1.reserve( nBestInliers );
		list2.reserve( nBestInliers );
		
		// iterator to the indices...
		typename std::vector< std::size_t >::const_iterator iti = iBestInliers.begin();
		const typename std::vector< std::size_t >::const_iterator itiEnd = iBestInliers.end();
		for ( ; iti < itiEnd; ++iti )
		{
			InputIterator1 it1 ( iBegin1 );
			InputIterator2 it2 ( iBegin2 );
			std::advance( it1, (*iti) );
			std::advance( it2, (*iti) );
			list1.push_back( *it1 );
			list2.push_back( *it2 );
		}
			

		typename RansacFunctor::Estimator()( result, list1.begin(), list1.end(), list2.begin(), list2.end() );
		OPT_LOG_DEBUG( "Estimated " << nBestInliers << " inlier after " << iRun + 1 << " iterations."  );
		return nBestInliers;
	}
	
	OPT_LOG_DEBUG( "RANSAC: Not enough inlier found" );
	return 0;
}

/**
 * RANSAC algorithm (for two-parameter problems)
 * @param result
 * @return 0 (failure) or number of inlier on success
 */
template< class Result, class Param1, class Param2, class Estimator, class Evaluator >
std::size_t ransac( Result& result, const std::vector< Param1 >& paramList1, const std::vector< Param2 >& paramList2,
	double fThreshold, std::size_t nSetSize, std::size_t nMinInlier, std::size_t nMinRuns, std::size_t nMaxRuns, 
	const Estimator& estimator, const Evaluator& evaluator, std::vector< bool >* pInliers = 0 )
{
	OPT_LOG_DEBUG( "RANSAC with " << paramList1.size() << " parameters, " << nMinInlier << " inlier required" );
	
	// set of inlier
	std::vector< bool > bInliers( paramList1.size() );

	// best so far
	std::size_t nBestInliers = 0;
	std::vector< bool > bBestInliers( paramList1.size() );
	
	std::size_t iRun ;
	for ( iRun = 0; iRun < nMaxRuns; iRun++ )
	{
		OPT_LOG_TRACE( "RANSAC iteration " << iRun + 1 );

		// try-catch-block will prevent crashes in case something becomes singular...
		try
		{
			// generate random set
			std::vector< Param1 > list1;
			std::vector< Param2 > list2;
			list1.reserve( nSetSize );
			list2.reserve( nSetSize );
			
			for ( std::size_t i = nSetSize; i; i-- )
			{
				std::size_t iSelected = rand() % paramList1.size();
				list1.push_back( paramList1[ iSelected ] );
				list2.push_back( paramList2[ iSelected ] );
				// TODO: check for double selections...
			}
			
			// compute hypothesis
			Result hypothesis;
			estimator( hypothesis, list1, list2 );
			
			// count inlier
			std::size_t nInlier = 0;
			double fInlierDist = 0;
			for ( std::size_t i = 0; i < paramList1.size() && int( paramList1.size() - i ) >= int( nMinInlier - nInlier ); i++ )
			{
				double d = evaluator( hypothesis, paramList1[ i ], paramList2[ i ] );
				if ( ( bInliers[ i ] = ( d < fThreshold ) ) == true )
				{
					nInlier++;
					fInlierDist += d;
				}
			}
			
			OPT_LOG_TRACE( nInlier << " inlier, avg dist=" << fInlierDist / nInlier );

			// save inlier set if we reached the required number or are better than a previous run
			if ( nInlier >= nMinInlier && nInlier > nBestInliers )
			{
				nBestInliers = nInlier;
				bBestInliers = bInliers;
			}

			// stop after nMinRun iterations if the required number of inlier was found
			if ( nBestInliers >= nMinInlier && iRun + 1 >= nMinRuns )
				break;
		}
#ifdef OPTIMIZATION_LOGGING
		catch ( const std::runtime_error& e )
		{ OPT_LOG_DEBUG( "RANSAC: caught exception: " << e.what() ); }
#else
		catch ( const std::runtime_error& )
		{}
#endif
	}

	if ( nBestInliers >= nMinInlier )
	{
		// compute final result
		std::vector< Param1 > list1;
		std::vector< Param2 > list2;
		list1.reserve( nBestInliers );
		list2.reserve( nBestInliers );
		for ( std::size_t i = 0; i < paramList1.size(); i++ )
			if ( bBestInliers[ i ] )
			{
				list1.push_back( paramList1[ i ] );
				list2.push_back( paramList2[ i ] );
			}
			
		estimator( result, list1, list2 );
		if ( pInliers )
			*pInliers = bBestInliers;
		OPT_LOG_DEBUG( iRun + 1 << " iterations, " << nBestInliers << " inlier" );

		return nBestInliers;
	}
	
	OPT_LOG_DEBUG( "RANSAC: Not enough inlier found" );
	return 0;
}

}}} // namespace Ubitrack::Math::Optimization

#endif // __UBITRACK_MATH_OPTIMIZATION_RANSAC_INCLUDED__

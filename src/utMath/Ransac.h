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
 */


#ifndef __Ubitrack_Math_Ransac_h_INCLUDED__
#define __Ubitrack_Math_Ransac_h_INCLUDED__

#include <utCore.h>
#include "Optimization.h"

namespace Ubitrack { namespace Math {

/**
 * RANSAC algorithm (for two-parameter problems)
 * @param result
 * @return 0 (failure) or number of inliers on success
 */
template< class Result, class Param1, class Param2, class Estimator, class Evaluator >
unsigned Ransac( Result& result, const std::vector< Param1 >& paramList1, const std::vector< Param2 >& paramList2,
	double fThreshold, unsigned nSetSize, unsigned nMinInliers, unsigned nMinRuns, unsigned nMaxRuns, 
	const Estimator& estimator, const Evaluator& evaluator, std::vector< bool >* pInliers = 0 )
{
	OPT_LOG_DEBUG( "RANSAC with " << paramList1.size() << " parameters, " << nMinInliers << " inliers required" );
	
	// set of inliers
	std::vector< bool > bInliers( paramList1.size() );

	// best so far
	unsigned nBestInliers = 0;
	std::vector< bool > bBestInliers( paramList1.size() );
	
	unsigned iRun ;
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
			
			for ( unsigned i = nSetSize; i; i-- )
			{
				unsigned iSelected = rand() % paramList1.size();
				list1.push_back( paramList1[ iSelected ] );
				list2.push_back( paramList2[ iSelected ] );
				// TODO: check for double selections...
			}
			
			// compute hypothesis
			Result hypothesis;
			estimator( hypothesis, list1, list2 );
			
			// count inliers
			unsigned nInliers = 0;
			double fInlierDist = 0;
			for ( unsigned i = 0; i < paramList1.size() && int( paramList1.size() - i ) >= int( nMinInliers - nInliers ); i++ )
			{
				double d = evaluator( hypothesis, paramList1[ i ], paramList2[ i ] );
				if ( ( bInliers[ i ] = ( d < fThreshold ) ) == true )
				{
					nInliers++;
					fInlierDist += d;
				}
			}
			
			OPT_LOG_TRACE( nInliers << " inliers, avg dist=" << fInlierDist / nInliers );

			// save inlier set if we reached the required number or are better than a previous run
			if ( nInliers >= nMinInliers && nInliers > nBestInliers )
			{
				nBestInliers = nInliers;
				bBestInliers = bInliers;
			}

			// stop after nMinRun interations if the required number of inliers was found
			if ( nBestInliers >= nMinInliers && iRun + 1 >= nMinRuns )
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

	if ( nBestInliers >= nMinInliers )
	{
		// compute final result
		std::vector< Param1 > list1;
		std::vector< Param2 > list2;
		list1.reserve( nBestInliers );
		list2.reserve( nBestInliers );
		for ( unsigned i = 0; i < paramList1.size(); i++ )
			if ( bBestInliers[ i ] )
			{
				list1.push_back( paramList1[ i ] );
				list2.push_back( paramList2[ i ] );
			}
			
		estimator( result, list1, list2 );
		if ( pInliers )
			*pInliers = bBestInliers;
		OPT_LOG_DEBUG( iRun + 1 << " iterations, " << nBestInliers << " inliers" );

		return nBestInliers;
	}
	
	OPT_LOG_DEBUG( "RANSAC: Not enough inliers found" );
	return 0;
}

} } // namespace Ubitrack::Math

#endif

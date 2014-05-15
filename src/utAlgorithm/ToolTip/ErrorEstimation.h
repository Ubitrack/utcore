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
 * @ingroup tracking_algorithms
 * @file
 * Error calculation for tooltip/hotspot calibration.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGROITHM_TOOLTIP_ERROR_ESTIMATION_H_INCLUDED__
#define __UBITRACK_ALGROITHM_TOOLTIP_ERROR_ESTIMATION_H_INCLUDED__

// Ubitrack
#include <utCore.h>
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Blas1.h>

// std
#include <utility> // std::pair
#include <numeric> // std::accumulate
#include <algorithm> // std::transform

namespace Ubitrack { namespace Algorithm { namespace ToolTip {


///  @internal functor object to calculate the resulting error
template< typename T >
struct ErrorFunction
{
protected:
	Math::Vector< T, 3 > tip;
	Math::Vector< T, 3 > offset;
	
public:
	ErrorFunction( const Math::Vector< T, 3 > & pw, const Math::Vector< T, 3 > & pm )
		: tip( pw )
		, offset( pm )
		{}
		
	T operator()( const Math::Pose& pose )
	{
		const Math::Vector< T, 3 > vec = tip - ( pose * offset );
		return Ubitrack::Math::norm_2( vec );
	}
	
	T operator()( const T sum, const Math::Pose& pose )
	{
		return sum + operator()( pose );
	}
};

/// @internal function to calculate the normal error of the tool tip calibration
template< typename T, typename InputIterator >
std::pair< T, T > estimatePosition3DError_6D( const Math::Vector< T, 3 >& pw
	, const InputIterator iBegin
	, const InputIterator iEnd
	, const Math::Vector< T, 3 >& pm )
{
	{	// old version
		// T err( 0 );
		// //err = std::accumulate( poses.begin(), poses.end(), err, ErrorFunction< T >( pw, pm ) );
		// err /= n;
	}

	// calculate distance errors in advance and reuse them later
	const std::size_t n ( std::distance(  iBegin, iEnd ) );
	std::vector< T > distance_error;
	distance_error.reserve( n );
	std::transform( iBegin, iEnd, std::back_inserter( distance_error ), ErrorFunction< T >( pw, pm ) );
	

	T err( 0 );
	{	// mean error
		err = std::accumulate( distance_error.begin(), distance_error.end(), err );
		err /= n;
	}

	
	T stdDev( 0 );
	{	// standard deviation
		for( typename std::vector< T >::const_iterator it ( distance_error.begin() ); it != distance_error.end(); ++it )
			stdDev += pow( *it - err, 2 );
		
		stdDev /= (n-1);
		stdDev = std::sqrt( stdDev );
	}
	return std::make_pair( err, stdDev ) ;
}
	
}}} // namespace Ubitrack::Algorithm::ToolTip


#endif //__UBITRACK_ALGROITHM_TOOLTIP_ERROR_ESTIMATION_H_INCLUDED__

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
 * Error calculation for absolute orientation problem
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_ERROR_ESTIMATION_H_INCLUDED__
#define __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_ERROR_ESTIMATION_H_INCLUDED__

// Ubitrack
#include <utCore.h>
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Blas1.h>
#include <utMath/Geometry/PointTransformation.h>

// std
#include <numeric> // std::accumulate
#include <algorithm> // std::transform

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation3D3D {


///  @internal functor object to calculate the resulting error
template< typename T >
struct ErrorFunction
{
protected:
	const Math::Pose& pose;
	
public:
	ErrorFunction( const Math::Pose &p )
		: pose( p )
		{}
		
	T operator()( const Math::Vector< T, 3 >& vecA, const Math::Vector< T, 3 >& vecB ) const
	{
		const Math::Vector< T, 3 > diffVector = vecA - ( pose * vecB );
		return Ubitrack::Math::norm_2( diffVector );
	}
};

/// @internal function to calculate the root-mean-square error of an absolute orientation
template< typename T, typename InputIterator >
T estimatePose6DResidual( const InputIterator itBegin1, const InputIterator itEnd1
		, const Math::Pose& pose, const InputIterator itBegin2, const InputIterator itEnd2 )
{
	// calculate distance errors in advance and reuse them later
	const std::size_t n ( std::distance(  itBegin1, itEnd1 ) );
	std::vector< T > distance_error;
	distance_error.reserve( n );
	std::transform( itBegin1, itEnd1, itBegin2, std::back_inserter( distance_error ), ErrorFunction< T >( pose ) );
	
	T err( 0 );
	{	// mean error, also residual, HartleyZissermann, p.136
		err = std::accumulate( distance_error.begin(), distance_error.end(), err );
		err /= (3*n);
		/// @todo check if a square root is necessary here:
		// err = std::sqrt( err );
	}
	return err;
}
	
}}} // namespace Ubitrack::Algorithm::PoseEstimation3D3D


#endif //__UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_ERROR_ESTIMATION_H_INCLUDED__

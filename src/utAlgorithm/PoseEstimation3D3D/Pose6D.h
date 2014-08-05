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
 * Implementation of Absolute Orientation (3D-3D-Pose-Estimation)
 *
 * @author Manuel Huber <huberma@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */

 
#ifndef __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_POSE_6D_INCLUDED__
#define __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_POSE_6D_INCLUDED__

#include "Rotation3D.h"
#include <utMath/Pose.h>
// #inlcude <utMath/Utils/PoseCast.h> <- will come in future..

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation3D3D {

/**
 * @internal
 * Calculate the absolute orientation problem.
 *
 * @param leftBegin iterator that points at the beginning of the 3D vectors in the left coordinate frame. Must be model of InputIterator and Math::Vector< double, 3 >*
 * @param leftEnd iterator that points at one after the last of the 3D vectors in the left coordinate frame. Must be model of InputIterator and Math::Vector< double, 3 >*
 * @param rightBegin iterator that points at the beginning of the 3D vectors in the right coordinate frame. Must be model of InputIterator and Math::Vector< double, 3 >*
 * @param rightEnd iterator that points at one after the last of the 3D vectors in the left coordinate frame. Must be model of InputIterator and Math::Vector< double, 3 >*
 * @return pose that describes the transformation of the left coordinate frame into the right coordinate frame.
 * @throws Util::Exception if lapack is not available, different number of samples for left and right side are given or the matrix N only has non-positive eigenvalues.
 */
template< typename InputIterator >
bool estimatePose6D_3D3D ( const InputIterator itBegin1, const InputIterator itEnd1
		, Math::Pose& pose, const InputIterator itBegin2, const InputIterator itEnd2
		, const typename std::iterator_traits< InputIterator >::value_type& leftCentroid
		, const typename std::iterator_traits< InputIterator >::value_type& rightCentroid )	
{
	/// @todo <-check that vector_type is a 3-vector
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type; 
	typedef typename vector_type::value_type value_type;

	{
		const std::size_t n = std::distance( itBegin1, itEnd1 );
		assert( n > 2u ); // <- 3D3D needs at least three pairwise distinct values
		const std::size_t n2 = std::distance( itBegin2, itEnd2 );
		assert( n == n2 ); // <- should contain equal amount of values
	}
	
	Math::Quaternion quat;
	if( !estimateRotation_3D3D ( itBegin1, itEnd1, quat, itBegin2, itEnd2, leftCentroid, rightCentroid ) )
		return false;
	
	const vector_type translation = leftCentroid - quat*rightCentroid;
	pose = Math::Pose ( quat, translation );
	
	return true;
}

/// @internal function overload that calculates the centroids as well
template< typename InputIterator >
bool estimatePose6D_3D3D ( const InputIterator itBegin1, const InputIterator itEnd1
		, Math::Pose& pose, const InputIterator itBegin2, const InputIterator itEnd2 )
{
	/// @todo <-check that vector_type is a 3-vector
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type; 

	const std::size_t n = std::distance( itBegin1, itEnd1 );

	const vector_type leftCentroid = std::accumulate( itBegin1, itEnd1, vector_type::zeros() ) / n;
	const vector_type rightCentroid = std::accumulate( itBegin2, itEnd2, vector_type::zeros() ) / n;

	return estimatePose6D_3D3D( itBegin1, itEnd1, pose, itBegin2, itEnd2, leftCentroid, rightCentroid );
}

} } } // namespace Ubitrack::Algorithm::PoseEstimation3D3D

#endif // __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_POSE_6D_INCLUDED__

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
 * A data structure that takes care of all relevant steps
 * of the hand eye calibration. The struct aims on easy usage
 * of the calibration procedure.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGORITHM_HANDEYE_STRUCT_H_INCLUDED__
#define __UBITRACK_ALGORITHM_HANDEYE_STRUCT_H_INCLUDED__

// Ubitrack
#include <utCore.h>
#include <utMath/Pose.h>
#include "DataSelection.h"
#include "DualQuaternion.h"

// std
#include <vector>
#include <iterator>

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation6D6D {


template< typename T >
struct HandEyeCalibration
{
	typedef T precision_type;
	typedef typename Math::Vector< precision_type, 8 > relative_pose_type;
	
	typedef typename std::vector< relative_pose_type > RelPoseListType;
	
	/// index pointing to the data, that was added latest
	std::size_t m_lastIndex;
	
	RelPoseListType relativePosesAllA; // the hands poses
	RelPoseListType relativePosesAllB; // the eye poses
	
	RelPoseListType relativePosesBestMatchA;
	RelPoseListType relativePosesBestMatchB;
	
	DataSelection< precision_type > m_DataSelection;
	
	HandEyeCalibration()
		: m_lastIndex( 0 )
		{}
		
	/// function that performs a hand eye calibration while it also updates internal structures that keep track of relative pose measurements and their best match to centroids	
	template< typename ResultPoseType, typename PoseTypeIn, template< typename Type, typename = std::allocator< Type > > class container_type >
	bool estimatePose6D( const container_type< PoseTypeIn >& posesEye, ResultPoseType& resultPose, const container_type< PoseTypeIn >& posesHand )
	{
		if( posesEye.size() < 3 )
			return false;
		
		if( posesEye.size() != posesHand.size() )
			return false;
		
		updateRelativePoses( posesEye, posesHand );
		resetBestMatch();
		return Algorithm::PoseEstimation6D6D::estimatePose6D_6D6D( relativePosesBestMatchB, resultPose, relativePosesBestMatchA );
	}
	
protected:	
	/// use this function to update to latest hand and eye relative pose measurements stored internally
	template< typename PoseType, template< typename Type, typename = std::allocator< Type > > class container_type >
	void updateRelativePoses( const container_type< PoseType >& eyePoseList, const container_type< PoseType >& handPoseList )
	{
		const std::size_t n = eyePoseList.size();
		assert( n == handPoseList.size() );
		updateRelativePoses< true >( handPoseList, relativePosesAllA );
		updateRelativePoses< false >( eyePoseList, relativePosesAllB );
		m_lastIndex = n;
	}
	
	/// @internal his function allocates new space for the relative pose measurements and triggers the next steps (== adding relative measurements)
	template< bool direction, typename PoseType, template< typename Type, typename = std::allocator< Type > > class container_type >
	void updateRelativePoses( const container_type< PoseType >& poseList, RelPoseListType& relPoseList ) const
	{
		// estimate the space needed for new pose measurements and reserve some space accordingly
		const std::size_t nRel = std::distance( relPoseList.begin(), relPoseList.end() );
		const std::size_t n = std::distance( poseList.begin(), poseList.end() );
		const std::size_t ni = n - m_lastIndex;
		const std::size_t nDiff = ni*m_lastIndex+(ni*(ni-1)/2);
		relPoseList.reserve( nRel + nDiff );
		
		updateRelativePoses< direction >( poseList.begin(), poseList.end(), m_lastIndex, std::back_inserter( relPoseList ) );
	}
	
	/// @internal this function updates the list of relative pose measurements from iterators and the last known position
	template< bool direction, typename InputIterator, typename OutputIterator >
	void updateRelativePoses( const InputIterator itBegin, const InputIterator itEnd, const std::size_t index, OutputIterator itOut ) const
	{
		// move pointer of incoming poses to newest measurement
		InputIterator itNew = itBegin;
		std::advance( itNew, index );
		
		for( ; itNew!=itEnd; ++itNew )
			for( InputIterator itOld = itBegin; itOld!=itNew; ++itOld )
				*itOut++ = ( Ubitrack::Math::relative_pose< typename RelPoseListType::value_type, direction >()( *itOld, *itNew ) );
	}
	
	/// @internal this function can be used to new relative pose measurements from a new pose measurement and a list of previous poses.
	template< bool direction, typename PoseType, template< typename Type, typename = std::allocator< Type > > class container_type >
	void updateRelativePoses( const container_type< PoseType >& oldPoseList, const PoseType& newPose, RelPoseListType& relPoseList ) const 
	{
		// allocate some more space to store the additional data
		const std::size_t nRel = std::distance( relPoseList.begin(), relPoseList.end() );
		const std::size_t nNew = std::distance( oldPoseList.begin(), oldPoseList.end() );
		relPoseList.reserve( nRel + nNew );
		//iterate over all previous poses to generate new relative poses
		typename container_type< PoseType >::const_iterator itBegin = oldPoseList.begin();
		const typename container_type< PoseType >::const_iterator itEnd = oldPoseList.end();
		for( ; itBegin != itEnd; ++itBegin )
			relPoseList.push_back( Ubitrack::Math::relative_pose< typename RelPoseListType::value_type, direction >()( *itBegin, newPose ) );
	}
	
	void resetBestMatchByAll()
	{
		relativePosesBestMatchA.clear();
		relativePosesBestMatchA.reserve( relativePosesAllA.size() );
		std::copy( relativePosesAllA.begin(), relativePosesAllA.end(), std::back_inserter( relativePosesBestMatchA ) );
		
		relativePosesBestMatchB.clear();
		relativePosesBestMatchB.reserve( relativePosesAllB.size() );
		std::copy( relativePosesAllB.begin(), relativePosesAllB.end(), std::back_inserter( relativePosesBestMatchB ) );
	}
	
	void resetBestMatch( const std::size_t m = 0 )
	{
		const std::size_t nRel = relativePosesAllA.size();
		assert( nRel == relativePosesAllB.size() );
		// std::cerr << "nRel " << nRel << "\n";
		// determine the number of k-means by simple rule (maybe not accurate)
		const std::size_t n = m ? m : static_cast< std::size_t >( sqrt( nRel/2. ) );
		
		// std::cerr << "n " << n << "\n";
		relativePosesBestMatchA.reserve( n );
		relativePosesBestMatchB.reserve( n );
		m_DataSelection.resetComparisonPoses( n, relativePosesAllA );
		m_DataSelection.getSelection( relativePosesAllA, std::back_inserter( relativePosesBestMatchA ) );
		m_DataSelection.getSelection( relativePosesAllB, std::back_inserter( relativePosesBestMatchB ) );
	}
};

}}} // namespace Ubitrack::Algorithm::PoseEstimation6D6D

#endif //__UBITRACK_ALGORITHM_HANDEYE_STRUCT_H_INCLUDED__

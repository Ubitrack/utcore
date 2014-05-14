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
 * Functions for ransac absolute orientation estimation .
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_RANSAC_H_INCLUDED__
#define __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_RANSAC_H_INCLUDED__


#include "Pose6D.h" // includes std::vector/Pose
#include <utCore.h>
#include <utMath/Blas1.h>
#include <utMath/Optimization/Ransac.h>

namespace Ubitrack { namespace Algorithm { namespace AbsoluteOrientation {

/**
 * @internal function object that provides estimation and evaluation
 * functions for a ransac absolute orientation estimation.
 */
template< typename T >
struct Ransac
{
public:

	typedef T value_type;
	
	/**
	 * @internal computes a 6d-pose from two given corresponding 3d-point sets via absolute orientation. 
	 */
	struct Estimator
	{
		public:
		
		template< typename InputIterator, typename ResultType >
		bool operator()( ResultType& resultPose, const InputIterator iBegin1, const InputIterator iEnd1, const InputIterator iBegin2, const InputIterator iEnd2 ) const
		{
			return estimatePose6D_3D3D( iBegin1, iEnd1, resultPose, iBegin2, iEnd2 );
		}
	};
	
	/**
	 * @internal computes euclidean distance of transformed point to original point
	 */
	struct Evaluator
	{
		public:
		
		template< typename PoseType, typename VectorType >
		T operator()( const PoseType &pose, const VectorType &vec1, const VectorType &vec2 ) const
		{
			const VectorType vecEstimated = pose * vec2;
			const VectorType diff( vec1 - vecEstimated );
			return Ubitrack::Math::norm_2( diff );
		}
	};
};

template< typename T, typename InputIterator, typename ResultType >
bool estimatePose6D_3D3D ( const InputIterator itBegin1, const InputIterator itEnd1
		, ResultType& pose
		, const InputIterator itBegin2, const InputIterator itEnd2
		, const Math::Optimization::RansacParameter< T >& params )
{
	const std::size_t inlier = Math::Optimization::ransac( itBegin1, itEnd1, itBegin2, itEnd2, pose, AbsoluteOrientation::Ransac< T >(), params  );
	return ( inlier > 0 );
}

UBITRACK_EXPORT bool estimatePose6D_3D3D( const std::vector< Math::Vector3f >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3f >& pointsB
	, const Math::Optimization::RansacParameter< float >& params );
	

UBITRACK_EXPORT bool estimatePose6D_3D3D( const std::vector< Math::Vector3d >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3d >& pointsB
	, const Math::Optimization::RansacParameter< double >& params );

}}} // namespace Ubitrack::Algorithm::AbsoluteOrientation

#endif //__UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_RANSAC_H_INCLUDED__
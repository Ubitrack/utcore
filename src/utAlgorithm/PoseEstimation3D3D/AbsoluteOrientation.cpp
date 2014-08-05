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
 * Implementation of  Absolute Orientation (3D-3D-Pose-Estimation)
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */

#include "AbsoluteOrientation.h"

#include "Scale.h"
#include "Rotation3D.h"
#include "Pose6D.h"
#include "Ransac.h"
#include "Optimization.h"

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation3D3D {

///  @internal old function call
Math::Pose calculateAbsoluteOrientation ( const std::vector< Math::Vector3d >& left,
										  const std::vector< Math::Vector3d >& right)
{
	Math::Pose pose;
	// attention flipped the input parameters for new order
	estimatePose6D_3D3D ( right.begin(), right.end(), pose, left.begin(), left.end() );
	return Math::Pose( pose );
}

double estimateScale_3D3D( const std::vector< Math::Vector3d >& m_left
	, const std::vector< Math::Vector3d >& m_right )
{
	return estimateScale_3D3D( m_left.begin(), m_left.end(), m_right.begin(), m_right.end() );
}

float estimateScale_3D3D( const std::vector< Math::Vector3f >& m_left
	, const std::vector< Math::Vector3f >& m_right )
{
	return estimateScale_3D3D( m_left.begin(), m_left.end(), m_right.begin(), m_right.end() );
}

bool estimatePose6D_3D3D( const std::vector< Math::Vector3d >& points3dA
	, Math::Pose& pose
	, const std::vector< Math::Vector3d >& points3dB )
{
	return estimatePose6D_3D3D( points3dA.begin(), points3dA.end(), pose, points3dB.begin(), points3dB.end() );
}

bool estimatePose6D_3D3D( const std::vector< Math::Vector3f >& points3dA
	, Math::Pose& pose, const std::vector< Math::Vector3f >& points3dB )
{
	return estimatePose6D_3D3D( points3dA.begin(), points3dA.end(), pose, points3dB.begin(), points3dB.end() );
}

bool estimatePose6D_3D3D( const std::vector< Math::Vector3f >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3f >& pointsB
	, const Math::Optimization::RansacParameter< float >& params )
{
	return estimatePose6D_3D3D( pointsA.begin(), pointsA.end(), pose, pointsB.begin(), pointsB.end(), params );
}

bool estimatePose6D_3D3D( const std::vector< Math::Vector3d >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3d >& pointsB
	, const Math::Optimization::RansacParameter< double >& params )
{
	return estimatePose6D_3D3D( pointsA.begin(), pointsA.end(), pose, pointsB.begin(), pointsB.end(), params );
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA
	, Math::Matrix3x3d& mat , const std::vector< Math::Vector3d >& points3dB )
{
	return estimateRotation_3D3D ( points3dA.begin(), points3dA.end(), mat, points3dB.begin(), points3dB.end() );
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3f >& points3dA
	, Math::Matrix3x3f& mat, const std::vector< Math::Vector3f >& points3dB )
{
	return estimateRotation_3D3D ( points3dA.begin(), points3dA.end(), mat, points3dB.begin(), points3dB.end() );
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA
	, Math::Quaternion& quat, const std::vector< Math::Vector3d >& points3dB )
{
	return estimateRotation_3D3D ( points3dA.begin(), points3dA.end(), quat, points3dB.begin(), points3dB.end() );
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3f >& points3dA
	, Math::Quaternion& quat, const std::vector< Math::Vector3f >& points3dB )
{
	return estimateRotation_3D3D ( points3dA.begin(), points3dA.end(), quat, points3dB.begin(), points3dB.end() );
}

/// @internal absolute orientation using non-linear optimization for \c float
bool estimatePose6D_3D3D( const std::vector< Math::Vector3f >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3f >& pointsB
	, const Math::Optimization::OptTerminate& criteria )
{
	return estimatePose6D_3D3D( pointsA.begin(), pointsA.end(), pose, pointsB.begin(), pointsB.end(), criteria );

}
/// @internal absolute orientation using non-linear optimization for \c double
bool estimatePose6D_3D3D( const std::vector< Math::Vector3d >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3d >& pointsB
	, const Math::Optimization::OptTerminate& criteria )
{
	return estimatePose6D_3D3D( pointsA.begin(), pointsA.end(), pose, pointsB.begin(), pointsB.end(), criteria );
}

} } } // namespace Ubitrack::Algorithm::PoseEstimation3D3D


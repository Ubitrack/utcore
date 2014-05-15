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
 * Implements all function calls for tooltip/hotspot calibration to generate
 * an object file providinig all available methods for standard datatypes.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */ 

// Ubitrack
#include "TipCalibration.h"
#include "ErrorEstimation.h"
#include "LeastSquare.h"
#include "Optimization.h"
#include "Ransac.h"


namespace Ubitrack { namespace Algorithm { namespace ToolTip {

/// @internal old function signature that performs the calibration
void tipCalibration( const std::vector< Math::Pose >& poses, 
	Math::Vector< double, 3 >& pm, Math::Vector< double, 3 >& pw )
{
	estimatePosition3D_6D< double >( pw, poses.begin(), poses.end(), pm );
}

/// @internal specialization of tooltip calibration for type \c float
bool estimatePosition3D_6D( Math::Vector< float, 3 >& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector< float, 3 >& pm )
{
	return estimatePosition3D_6D< float >( pw, poses.begin(), poses.end(), pm );
}

/// @internal specialization of tooltip calibration for type \c double
bool estimatePosition3D_6D( Math::Vector< double, 3 >& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector< double, 3 >& pm )
{
	return estimatePosition3D_6D< double >( pw, poses.begin(), poses.end(), pm );
}

/// @internal specialization of tooltip error estimation for type \c float
std::pair< float, float > estimatePosition3DError_6D( const Math::Vector3f& pw
	, const std::vector< Math::Pose >& poses 
	, const Math::Vector3f& pm )
{
	return estimatePosition3DError_6D( pw, poses.begin(), poses.end(), pm );
}

/// @internal specialization of tooltip error estimation for type \c double
std::pair< double, double > estimatePosition3DError_6D( const Math::Vector3d& pw
	, const std::vector< Math::Pose >& poses 
	, const Math::Vector3d& pm )
{
	return estimatePosition3DError_6D( pw, poses.begin(), poses.end(), pm );
}

/// @internal specialization of non-linearly optimized tooltip calibration for type \c float
bool estimatePosition3D_6D( Math::Vector3f& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3f& pm
	, const Math::Optimization::OptTerminate& criteria )
{
	return estimatePosition3D_6D( pw, poses.begin(), poses.end(), pm, criteria );
}

/// @internal specialization of non-linearly optimized tooltip calibration for type \c double
bool estimatePosition3D_6D( Math::Vector3d& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3d& pm
	, const Math::Optimization::OptTerminate& criteria )
{
	return estimatePosition3D_6D( pw, poses.begin(), poses.end(), pm, criteria );
}


/// @internal specialization of tooltip ransac calibration for type \c float
bool estimatePosition3D_6D( Math::Vector3f& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3f& pm
	, const Math::Optimization::RansacParameter< float >& params )
{
	return estimatePosition3D_6D( pw, poses.begin(), poses.end(), pm, params );
}
	
/// @internal specialization of tooltip ransac calibration for type \c double
bool estimatePosition3D_6D( Math::Vector3d& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3d& pm
	, const Math::Optimization::RansacParameter< double >& params )
{
	return estimatePosition3D_6D( pw, poses.begin(), poses.end(), pm, params );
}

}}} // namespace Ubitrack::Algorithm::ToolTip

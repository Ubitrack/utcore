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
 * Functions for tooltip/hotspot calibration.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */ 

#ifndef __UBITRACK_ALGROITHM_TOOLTIP_CALIBRATION_H_INCLUDED__
#define __UBITRACK_ALGROITHM_TOOLTIP_CALIBRATION_H_INCLUDED__

// Ubitrack
#include <utCore.h>
#include <utMath/Pose.h>
#include <utMath/Vector.h> //incldues std::vector

// std
#include <utility> // std::pair

namespace Ubitrack { namespace Algorithm { namespace ToolTip {

/**
 * @ingroup tracking_algorithms
 * Computes the tooltip/hotspot calibration.
 *
 * The routine solves the following equation system using a least-square solution
 * , given a list of body @f$ i @f$ poses (R_i, t_i):
 *  @f$ (R_i -I) (p_m p_w) = -t_i @f$
 *
 * A description of the algorithm was published by Tuceryan et al. in their article
 * "Calibration requirements and procedures for a monitor-based augmented reality system"
 * in 1995 ( @cite tuceryan1995calibration ).
 *
 * @verbatim
@article{tuceryan1995calibration,
  title={Calibration requirements and procedures for a monitor-based augmented reality system},
  author={Tuceryan, Mihran and Greer, Douglas S. and Whitaker, Ross T. and Breen, David E. and Crampton, Chris and Rose, Eric and Ahlers, Klaus H},
  journal={Visualization and Computer Graphics, IEEE Transactions on},
  volume={1},
  number={3},
  pages={255--273},
  year={1995},
  publisher={IEEE}
} @endverbatim
 *
 *
 * @param pm returns the constant point in body coordinates
 * @param poses the list of poses that
 * @param pw returns the constant point in world coordinates
 */	
UBITRACK_EXPORT bool estimatePosition3D_6D( Math::Vector3f& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3f& pm );
	
UBITRACK_EXPORT bool estimatePosition3D_6D( Math::Vector3d& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3d& pm );
	
/**
 * @ingroup tracking_algorithms
 * Computes the gaussian error of a tip/hotspot calibration.
 *
 * The routine calculates the mean error as
 * @f$ \frac{\sum_i=0^n=N | p_w - R_i p_m + t_i |_2}{n} @f$ 
 *
 * @param pw the constant point(tooltip/hotspot) in world coordinates
 * @param poses the container of poses describing the orbit around the (tooltip/hotspot)
 * @param pm the constant point(tooltip/hotspot) in body coordinates
 * @return returns a pair of values including the mean value and the standard deviation of the error
 */
UBITRACK_EXPORT std::pair< float, float > estimatePosition3DError_6D( const Math::Vector3f& pw
	, const std::vector< Math::Pose >& poses 
	, const Math::Vector3f& pm );

UBITRACK_EXPORT std::pair< double, double > estimatePosition3DError_6D( const Math::Vector3d& pw
	, const std::vector< Math::Pose >& poses 
	, const Math::Vector3d& pm );


/// old version of function call, please use instead: \c bool estimatePosition3D_6D( Math::Vector3f& pw, const std::vector< Math::Pose >& poses , Math::Vector3f& pm );
UBITRACK_EXPORT void tipCalibration( const std::vector< Math::Pose >& poses, 
	Math::Vector3d& pm, Math::Vector3d& pw );
	
}}} // namespace Ubitrack::Algorithm::ToolTip


#endif //__UBITRACK_ALGROITHM_TOOLTIP_CALIBRATION_H_INCLUDED__

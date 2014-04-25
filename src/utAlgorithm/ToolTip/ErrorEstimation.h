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
#include <utMath/Vector.h> //includes std::vector

// std
#include <utility> // std::pair

namespace Ubitrack { namespace Algorithm { namespace ToolTip {

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

	
}}} // namespace Ubitrack::Algorithm::ToolTip


#endif //__UBITRACK_ALGROITHM_TOOLTIP_ERROR_ESTIMATION_H_INCLUDED__

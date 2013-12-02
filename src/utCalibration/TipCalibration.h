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
 * Functions for tip/hotspot calibration.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_CALIBRATION_TIPCALIBRATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_TIPCALIBRATION_H_INCLUDED__



#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/Pose.h>
#include <vector>

namespace Ubitrack { namespace Calibration {

/**
 * @ingroup tracking_algorithms
 * Computes the tip/hotspot calibration.
 *
 * The routine solves the following equation system, given a list of body poses (R_i, t_i):
 *  (R_i -I) (p_m p_w) = -t_i
 *
 * @tparam T type of parameters ( e.g. \c double or \c float )
 * @param poses the list of poses
 * @param pm returns the constant point in body coordinates
 * @param pw returns the constant point in world coordinates
 */
 void tipCalibrationImpl( const std::vector< Math::Pose >& poses, 
	Math::Vector< 3, double >& pm, Math::Vector< 3, double >& pw );

/**
 * @ingroup tracking_algorithms
 * Computes the tip/hotspot calibration.
 *
 * The routine solves the following equation system, given a list of body poses (R_i, t_i):
 *  (R_i -I) (p_m p_w) = -t_i
 *
 * @param poses the list of poses
 * @param pm returns the constant point in body coordinates
 * @param pw returns the constant point in world coordinates
 */
UBITRACK_EXPORT void tipCalibration( const std::vector< Math::Pose >& poses, 
	Math::Vector< 3, double >& pm, Math::Vector< 3, double >& pw );

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

#endif

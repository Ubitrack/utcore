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
 * Functions for HandEye Calibration 
 *
 * @author Daniel Muhra <muhra@in.tum.de>
 */ 

#ifndef __UBITRACK_CALIBRATION_HANDEYE_CALIBRATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_HANDEYE_CALIBRATION_H_INCLUDED__



#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/Matrix.h>
#include <utMath/Pose.h>
#include <vector>

namespace Ubitrack { namespace Calibration {

/**
 * @ingroup tracking_algorithms
 * Computes a Pose using a Hand-Eye-Calibration.
 * This Implementation uses the algorithm described by Tsai - Lenz.
 * You can use either a series of double or float matrices.
 * Using of pose instead of matrices is also possible.
 * Eye and hand vectors must contain the same number of objects
 *
 * @param hand vector containing a series of hand (marker) poses in the global (tracker) coordinate system
 * @param eye vector containing a series of eye (camera) poses in the eye coordinate system
 * @return the transformation between eye and hand (e.g. camera and attached marker)
 */
UBITRACK_EXPORT Math::Pose performHandEyeCalibration ( const std::vector< Math::Matrix< float, 4, 4 > >& hand,  const std::vector< Math::Matrix< float, 4, 4 > >& eye, bool bUseAllPairs = true );

UBITRACK_EXPORT Math::Pose performHandEyeCalibration ( const std::vector< Math::Matrix< double, 4, 4 > >& hand,  const std::vector< Math::Matrix< double, 4, 4 > >& eye, bool bUseAllPairs = true );

UBITRACK_EXPORT Math::Pose performHandEyeCalibration ( const std::vector< Math::Pose >& hand,  const std::vector< Math::Pose >& eye, bool bUseAllPairs = true );


}} // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

#endif

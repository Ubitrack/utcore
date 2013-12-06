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
 * Functions for 2D-3D pose estimation.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_2D3DPOSEESTIMATION_HAGER_H_INCLUDED__
#define __UBITRACK_CALIBRATION_2D3DPOSEESTIMATION_HAGER_H_INCLUDED__


#include <vector>
#include <utCore.h>
#include <utMath/Matrix.h>
#include <utMath/Pose.h>

namespace Ubitrack { namespace Calibration {

#ifdef HAVE_LAPACK

/**
 * @ingroup tracking_algorithms
 * Hager's fast and globally convergent pose estimation @cite lu2000fast.
 * 
 *
 * Note: Also exists with \c double parameters.
 *
 * @param p the initial pose
 * @param p2D points in image coordinates
 * @param p3D points in object coordinates
 * @param cam camera intrinsics matrix
 * @param nIterations number of iterations
 * @return residual of the optimization
 */

UBITRACK_EXPORT bool estimatePose( Math::Pose& p,
	std::vector< Math::Vector< float, 3 > > p2D,
	const std::vector< Math::Vector< float, 3 > > p3D,
	const Math::Matrix< 3, 3, float >& cam,
	unsigned &nIterations,
	float &error );

UBITRACK_EXPORT bool estimatePose( Math::Pose& p,
	std::vector< Math::Vector< double, 3 > > p2D,
	const std::vector< Math::Vector< double, 3 > > p3D,
	const Math::Matrix< 3, 3, double >& cam ,
	unsigned &nIterations ,
	double &error );
	
#endif // HAVE_LAPACK
	
} } // namespace Ubitrack::Calibration

#endif

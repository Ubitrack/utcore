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

#include <utCore.h>
#include <utMath/Pose.h>

#include <vector>

namespace Ubitrack { namespace Calibration {

#ifdef HAVE_LAPACK

/**
 * @ingroup calibration tracking_algorithms
 * @brief An algorithm to determine a solution to the classic \b 2D-3D 
 * \b pose \b estimation problem in monocular vision scenarios from given \b 3D points
 * and corresponding \b 2D observations.
 *
 * This algorithm estimates a \b pose from given \b 2D and \b 3D point correspondences.
 * This problem is of interest since the beginning of photogrammetry and many
 * solution can be found that provide a solution to this problem.
 * The hereby implemented approach is based on
 * "Fast and Globally Convergent Pose Estimation from Video Images"
 * from Lu et al. in 2000 ( @cite lu2000fast ):
 *
 * @verbatim
@article{lu2000fast,
 title={Fast and globally convergent pose estimation from video images},
 author={Lu, Chien-Ping and Hager, Gregory D. and Mjolsness, Eric},
 journal={Pattern Analysis and Machine Intelligence, IEEE Transactions on},
 volume={22},
 number={6},
 pages={610--622},
 year={2000},
 publisher={IEEE} 
} @endverbatim
 *
 * 
 * Example use case:\n
 * @code
 * std::vector< Vector3d > points3d; // <- should be filled with 3D object points
 * std::vector< Vector2d > points2d; // <- should be filled with projected 3D points 
 * Math::Pose pose; // <- will be filled with values
 * estimatePose6D_2D3D( points2d, pose, points3d, 50, 1e-06 );
 * @endcode
 *
 * @attention : There is a version of this function ovlerloaded with \c float instead of \c double parameters.
 *
 * @param p2D points in (normalized) image coordinates (mabye you need to apply @f$ K^{-1} @f$ ) to the 2d points first
 * @param pose the (initial) pose and final result
 * @param p3D points in object coordinates
 * @param max_iter maximum number of allowed iterations
 * @param min_error the minimum change in error allowed to converge
 * @return flag that signs if the algorithm converged due to minimal error
 */
UBITRACK_EXPORT bool estimatePose6D_2D3D(  const std::vector< Math::Vector2d >& p2D, Math::Pose& pose,
	const std::vector< Math::Vector3d >& p3D, std::size_t &max_iter, double &min_error );

/// @internal overloaded function with float values
UBITRACK_EXPORT bool estimatePose6D_2D3D( const std::vector< Math::Vector2f >& p2D, Math::Pose& pose,
	const std::vector< Math::Vector3f >& p3D, std::size_t &nIterations, float &error );

#endif // HAVE_LAPACK
	
} } // namespace Ubitrack::Calibration

#endif

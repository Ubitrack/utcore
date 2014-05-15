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

// #include <utMeasurement/Measurement.h>

#include "../utCore.h"

#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Pose.h>

namespace Ubitrack { namespace Algorithm {

#ifdef HAVE_LAPACK

/**
 * @ingroup tracking_algorithms
 * @brief Performs a (classic) bundle adjustment
 *
 * This method performs an optimization on given (noisy) 3D point data and 
 * (noisy) camera poses from given 2D image observations of the 3D points
 * using a common bundle adjustment approach. 
 *
 * The function to be minimized is stated as \@f  min \sum_{ij} d( \hat{P}^{i} \hat{X}_j, x_{j}^{i})^{2} \@f 
 * The Jacobian in here is build up similar to figure A.1
 *  on http://www.cs.unc.edu/~marc/tutorial/node163.html#sec:subbundle .
 *
 * @param pts2D \c std::vector of observations, for each camera a new \c std::vector
 * @return camPoses \c std::vector of poses, defining the initial extrinsic camera orientations
 * @return pts3D \c std::vector of initial 3D points , basis of the 2D observations in 1st parameter
 */

UBITRACK_EXPORT void simpleBundleAdjustment( const std::vector< std::vector< Math::Vector2d > >& pts2D, std::vector< Math::Pose >& camPoses, std::vector< Math::Vector3d >& pts3D );

UBITRACK_EXPORT void simpleBundleAdjustment( const std::vector< std::vector< Math::Vector2f > >& pts2D, std::vector< Math::Pose >& camPoses, std::vector< Math::Vector3f >& pts3D);
	
#endif // HAVE_LAPACK

}} // namespace Ubitrack::Algorithm

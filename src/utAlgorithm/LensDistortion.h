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
 * Functions for lens distortion.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_CALIBRATION_LENSDISTORTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_LENSDISTORTION_H_INCLUDED__

#include <utCore.h>
#include <utMath/Matrix.h>
#include <utMath/Vector.h>

namespace Ubitrack { namespace Algorithm {

/**
 * projects a 3D point in camera coordinates onto the image plane, applying radial distortion.
 * @param p the 3D point in camera coordinates
 * @param dist distortion parameters
 * @param K 3x3 intrinsic matrix
 * @return projected and distorted point in image coordinates
 */
UBITRACK_EXPORT Math::Vector< float, 2 > projectWithDistortion( const Math::Vector< float, 3 >& p, const Math::Vector< float, 4 >& dist,
	const Math::Matrix< float, 3, 3 >& K );

UBITRACK_EXPORT Math::Vector< double, 2 > projectWithDistortion( const Math::Vector< double, 3 >& p, const Math::Vector< double, 4 >& dist,
	const Math::Matrix< double, 3, 3 >& K );


/**
 * apply lens distortion to a point already in image coordinates.
 * @param p 2D image point
 * @param dist distortion parameters
 * @param K 3x3 intrinsic matrix
 * @return distorted point in image coordinates
 */
UBITRACK_EXPORT Math::Vector< float, 2 > lensDistort( const Math::Vector< float, 2 >& p, const Math::Vector< float, 4 >& dist,
	const Math::Matrix< float, 3, 3 >& K );

UBITRACK_EXPORT Math::Vector< double, 2 > lensDistort( const Math::Vector< double, 2 >& p, const Math::Vector< double, 4 >& dist,
	const Math::Matrix< double, 3, 3 >& K );

#ifdef HAVE_LAPACK
	
/**
 * remove lens distortion to a point already in image coordinates.
 * @param p 2D distorted image point
 * @param dist distortion parameters
 * @param K 3x3 intrinsic matrix
 * @return undistorted point in image coordinates
 */
UBITRACK_EXPORT Math::Vector< float, 2 > lensUnDistort( const Math::Vector< float, 2 >& p, const Math::Vector< float, 4 >& dist,
	const Math::Matrix< float, 3, 3 >& K );

UBITRACK_EXPORT Math::Vector< double, 2 > lensUnDistort( const Math::Vector< double, 2 >& p, const Math::Vector< double, 4 >& dist,
	const Math::Matrix< double, 3, 3 >& K );

#endif
	
} } // namespace Ubitrack::Algorithm

#endif

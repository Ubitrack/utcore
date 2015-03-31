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
 * Functions for lens un/-distortion.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGORITHM_CAMERALENS_CORRECTION_H_INCLUDED__
#define __UBITRACK_ALGORITHM_CAMERALENS_CORRECTION_H_INCLUDED__

#include <vector>	// std::vector

#include <utCore.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/CameraIntrinsics.h>

namespace Ubitrack { namespace Algorithm { namespace CameraLens {

/**
 * apply lens distortion to a point given in image coordinates
 *
 * Distorts a 2d point radially and tangentially. The distortion is described by a two vectors
 * in addition the point is unprojected and projected from image coordinates to sensor coordinates and back
 * using the intrinsics matrix from the calibration parameters.
 *
 * a 6-vector containing the coefficients for the radial distortion ( k1, k2[, k3[, k4, k5, k6]]
 * a 2-vector containing the coefficients for the tangential distortion ( p1, p2 )
 * 
 * Formula for 2 radial distortion parameters
 *    x' = x + x( k1 * r^2 + k2 * r^4 ) + ( 2 * p1 * x * y + p2 * ( r^2 + 2 * x^2 ) )
 *    y' = y + y( k1 * r^2 + k2 * r^4 ) + ( 2 * p2 * x * y + p1 * ( r^2 + 2 * y^2 ) )
 *
 * Formula for 3 radial distortion parameters
 *    x' = x + x( k1 * r^2 + k2 * r^4 + k3 * r^6 ) + ( 2 * p1 * x * y + p2 * ( r^2 + 2 * x^2 ) )
 *    y' = y + y( k1 * r^2 + k2 * r^4 + k3 * r^6 ) + ( 2 * p2 * x * y + p1 * ( r^2 + 2 * y^2 ) )
 *
 * Formula for 6 radial distortion parameters
 *    x' = x( ( 1 + k1 * r^2 + k2 * r^4 + k3 * r^6 ) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6) ) + ( 2 * p1 * x * y + p2 * ( r^2 + 2 * x^2 ) )
 *    y' = x( ( 1 + k1 * r^2 + k2 * r^4 + k3 * r^6 ) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6) ) + ( 2 * p2 * x * y + p1 * ( r^2 + 2 * y^2 ) )
 * where
 *    r^2 = x^2 + y^2
 *
 * @attention : There are overloaded versions of this function for \c double precision types and vectors of 2d-images points
 *
 * @param intrinsics camera intrinsics parameters including 3x3 intrinsic matrix and distortion parameters
 * @param undistorted undistorted point in image coordinates
 * @param distorted distorted 2d image point
 */
UBITRACK_EXPORT void distort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& undistorted, Math::Vector2f& distorted );

/** 
 * @brief overloaded function \c distort with \c float parameters.
 *
 * For further information on this algorithm see void distort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& distorted, Math::Vector2f& undistorted );
 */
UBITRACK_EXPORT void distort( const Math::CameraIntrinsics< double >& intrinsics, const Math::Vector2d& undistorted, Math::Vector2d& distorted );

/** 
 * @brief overloaded function \c distort with \c for vectors of 2d points with single precision parameters.
 *
 * For further information on this algorithm see void distort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& distorted, Math::Vector2f& undistorted );
 */
UBITRACK_EXPORT void distort( const Math::CameraIntrinsics< float >& intrinsics, const std::vector< Math::Vector2f >& undistorted, std::vector< Math::Vector2f >& distorted );

/** 
 * @brief overloaded function \c distort with \c for vectors of 2d points with double precision parameters.
 *
 * For further information on this algorithm see void distort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& distorted, Math::Vector2f& undistorted );
 */
UBITRACK_EXPORT void distort( const Math::CameraIntrinsics< double >& intrinsics, const std::vector< Math::Vector2d >& undistorted, std::vector< Math::Vector2d >& distorted );


#ifdef HAVE_LAPACK

/**
 * remove lens distortion to a point 
 *
 * the point should be in image coordinates (pixels)
 *
 * Unistorts a 2d point radially and tangentially. The distortion is described by a two vectors
 * in addition the point is unprojected and projected from image coordinates to sensor coordinates and back
 * using the intrinsics matrix from the calibration parameters.
 *
 * a 6-vector containing the coefficients for the radial distortion ( k1, k2[, k3[, k4, k5, k6]]
 * a 2-vector containing the coefficients for the tangential distortion ( p1, p2 )
 * 
 * applies a non-linear optimization to the image points to approximate the undistorted image points. The optimization is based on the
 * distortion function explained in void distort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& undistorted, Math::Vector2f& distorted );
 *
 * @attention : There are overloaded versions of this function for \c double precision types and vectors of 2d-images points
 *
 * @param intrinsics camera intrinsics parameters including 3x3 intrinsic matrix and distortion parameters
 * @param distorted distorted 2d image point
 * @param undistorted undistorted point in image coordinates
 */
UBITRACK_EXPORT void undistort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& distorted, Math::Vector2f& undistorted );


/** 
 * @brief overloaded function \c undistort with \c float parameters.
 *
 * For further information on this algorithm see void undistort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& distorted, Math::Vector2f& undistorted );
 */
UBITRACK_EXPORT void undistort( const Math::CameraIntrinsics< double >& intrinsics, const Math::Vector2d& distorted, Math::Vector2d& undistorted );

/** 
 * @brief overloaded function \c undistort with \c for vectors of 2d points with single precision parameters.
 *
 * For further information on this algorithm see void undistort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& distorted, Math::Vector2f& undistorted );
 */
UBITRACK_EXPORT void undistort( const Math::CameraIntrinsics< float >& intrinsics, const std::vector< Math::Vector2f >& distorted, std::vector< Math::Vector2f >& undistorted );

/** 
 * @brief overloaded function \c undistort with \c for vectors of 2d points with double precision parameters.
 *
 * For further information on this algorithm see void undistort( const Math::CameraIntrinsics< float >& intrinsics, const Math::Vector2f& distorted, Math::Vector2f& undistorted );
 */
UBITRACK_EXPORT void undistort( const Math::CameraIntrinsics< double >& intrinsics, const std::vector< Math::Vector2d >& distorted, std::vector< Math::Vector2d >& undistorted );
	
#endif
	
}}} // namespace Ubitrack::Algorithm::CameraLens

#endif	//__UBITRACK_ALGORITHM_CAMERALENS_CORRECTION_H_INCLUDED__

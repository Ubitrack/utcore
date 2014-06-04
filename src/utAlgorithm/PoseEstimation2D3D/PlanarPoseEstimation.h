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
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_ALGORITHM_PLANAR_2D3D_POSE_ESTIMATION_H_INCLUDED__
#define __UBITRACK_ALGORITHM_PLANAR_2D3D_POSE_ESTIMATION_H_INCLUDED__



#include <utCore.h>
#include <utMath/Matrix.h>
#include <utMath/Pose.h>
#include <utMath/ErrorPose.h>


namespace Ubitrack { namespace Algorithm { namespace PoseEstimation2D3D {
/**
 * @ingroup tracking_algorithms
 * Computes a pose given a homography.
 * Note: Also exists with \c double parameters.
 *
 * @param H the homography that maps world to image coordinates
 * @param invK the inverse of the camera intrinsics matrix K that maps camera coordinates to image coordinates
 * @return the camera pose
 */
UBITRACK_EXPORT Math::Pose poseFromHomography( const Math::Matrix< float, 3, 3 >& H, const Math::Matrix< float, 3, 3 >& invK );

UBITRACK_EXPORT Math::Pose poseFromHomography( const Math::Matrix< double, 3, 3 >& H, const Math::Matrix< double, 3, 3 >& invK );


#ifdef HAVE_LAPACK

/**
 * @ingroup tracking_algorithms
 * Optimize a pose with levenberg-marquardt given 2D-3D correspondences.
 * Note: Also exists with \c double parameters.
 *
 * @param p the initial pose
 * @param p2D points in image coordinates
 * @param p3D points in object coordinates
 * @param cam camera intrinsics matrix
 * @param nIterations number of levenberg-marquardt iterations
 * @return residual of the optimization
 */
UBITRACK_EXPORT float optimizePose( Math::Pose& p, const std::vector< Math::Vector< float, 2 > >& p2D, 
	const std::vector< Math::Vector< float, 3 > >& p3D, const Math::Matrix< float, 3, 3 >& cam, 
	const std::size_t nIterations = 6 );

UBITRACK_EXPORT double optimizePose( Math::Pose& p, const std::vector< Math::Vector< double, 2 > >& p2D, 
	const std::vector< Math::Vector< double, 3 > >& p3D, const Math::Matrix< double, 3, 3 >& cam,
	const std::size_t nIterations = 6 );

	
/**
 * @ingroup tracking_algorithms
 * Computes the covariance of a pose created from observations of known 3D points by a single camera.
 * For an explanation of the covariance matrix, see \c ErrorPose.
 * Note: Also exists with \c double parameters.
 *
 * @param p the already computed pose whose error is to be estimated (extrinsic camera parameters)
 * @param p3D the points in object coordinates
 * @param cam camera intrinsics matrix
 * @param imageError variance of the measurement error in the image plane
 */
UBITRACK_EXPORT Math::Matrix< float, 6, 6 > singleCameraPoseError( const Math::Pose& p, const std::vector< Math::Vector< float, 3 > >& p3D, 
	const Math::Matrix< float, 3, 3 >& cam, float imageError );
	
UBITRACK_EXPORT Math::Matrix< double, 6, 6 > singleCameraPoseError( const Math::Pose& p, const std::vector< Math::Vector< double, 3 > >& p3D, 
	const Math::Matrix< double, 3, 3 >& cam, double imageError );

/**
 * @ingroup tracking_algorithms
 * Computes the covariance of a pose created from observations of known 3D points by a multiple cameras, 
 * e.g. the ART DTrack tracking system.
 * For an explanation of the covariance matrix, see \c ErrorPose.
 * Note: Also exists with \c double parameters.
 *
 * @param p the already computed pose whose error is to be estimated
 * @param p3D the points in object coordinates
 * @param cameras 3x4 projection matrices of the cameras containing both intrinsic and extrinsic parameters
 * @param observations a list of tuples (i_p, i_c) meaning that camera i_c has seen point i_p, where
 *   i_p and i_c are the indeces in the \c p3D and \c cameras arrays.
 * @param imageError variance of the measurement error in the image plane, assumed to be the same for all
 *   cameras and image points.
 */
UBITRACK_EXPORT Math::Matrix< float, 6, 6 > multipleCameraPoseError( const Math::Pose& p, 
	const std::vector< Math::Vector< float, 3 > >& p3D, 
	const std::vector< Math::Matrix< float, 3, 4 > >& cameras, 
	const std::vector< std::pair< std::size_t, std::size_t > > observations, 
	float imageError );
	
UBITRACK_EXPORT Math::Matrix< double, 6, 6 > multipleCameraPoseError( const Math::Pose& p, 
	const std::vector< Math::Vector< double, 3 > >& p3D, 
	const std::vector< Math::Matrix< double, 3, 4 > >& cameras, 
	const std::vector< std::pair< std::size_t, std::size_t > > observations, 
	double imageError );

	
/**
 * Initialization type needed for computePose() method.
 * Use \c NONPLANAR_PROJECTION only in case you are sure that points are not coplanar.
 */
UBITRACK_EXPORT typedef enum InitializationMethod {
	PLANAR_HOMOGRAPHY,
	NONPLANAR_PROJECTION
} InitializationMethod_t;


/**
 * @ingroup tracking_algorithms
 * Computes a pose given 2D-3D point correspondences
 * @param p2D points in image coordinates
 * @param p3D points in object coordinates
 * @param cam camera intrinsics matrix
 * @param initMethod Method used for initialization of the non-linear optimization
 */
UBITRACK_EXPORT Math::ErrorPose computePose( 
		const std::vector< Math::Vector< double, 2 > >& p2d,
		const std::vector< Math::Vector< double, 3 > >& p3d,
		const Math::Matrix< double, 3, 3 >& cam,
		bool optimize = true,
		enum InitializationMethod initMethod = (enum InitializationMethod)PLANAR_HOMOGRAPHY
	);
	
/**
 * @ingroup tracking_algorithms
 * Computes a pose given 2D-3D point correspondences
 * @param p2D points in image coordinates
 * @param p3D points in object coordinates
 * @param cam camera intrinsics matrix
 * @param residual reprojection error im image coordinates
 * @param initMethod Method used for initialization of the non-linear optimization
 */
UBITRACK_EXPORT Math::ErrorPose computePose( 
		const std::vector< Math::Vector< double, 2 > >& p2d,
		const std::vector< Math::Vector< double, 3 > >& p3d,
		const Math::Matrix< double, 3, 3 >& cam,
		double& residual,
		bool optimize = true,
		enum InitializationMethod initMethod = (enum InitializationMethod)PLANAR_HOMOGRAPHY		
	);
	
#endif // HAVE_LAPACK
	
} } } // namespace Ubitrack::Algorithm::PoseEstimation2D3D

#endif // __UBITRACK_ALGORITHM_PLANAR_2D3D_POSE_ESTIMATION_H_INCLUDED__

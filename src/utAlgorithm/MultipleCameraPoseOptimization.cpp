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
 * @ingroup calibration
 * @file
 * 2D-3D pose optimization component for multiple-camera systems.
 * Requires an initial pose!
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
 
#include "MultipleCameraPoseOptimization.h"



//#define OPTIMIZATION_LOGGING

// get a logger
#include <log4cpp/Category.hh>
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Calibration.2D6DPoseEstimation" ) );
//static log4cpp::Category& optLogger( log4cpp::Category::getInstance( "Ubitrack.Calibration.2D6DPoseEstimation.LM" ) );


#include <utMath/Optimization/LevenbergMarquardt.h>
#include "PoseEstimation2D3D/2D3DPoseEstimation.h"
#include <utUtil/Exception.h>



namespace Ubitrack { namespace Algorithm {

// In future, this file should not be compiled at all if lapack is not available
#ifdef HAVE_LAPACK

std::pair < Math::ErrorPose , double > 
	multipleCameraEstimatePose (
	const std::vector < Math::Vector< double, 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector< double, 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar< double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< double, 3, 3 > >& camMatrices,
	const int minCorrespondences,
	bool hasInitialPoseProvided,
	Math::Pose initialPose,
	int startIndex,
	int endIndex)
{
	if (endIndex == -1)
		endIndex = static_cast<int>( points3d.size() ) - 1 ;

	namespace ublas = boost::numeric::ublas;
	const std::size_t numberCameras ( points2dWeights.size() );


	// observation vector list which is used by the objective function
	std::vector< std::pair< std::size_t, std::size_t > > observations;

	// local 3dpoints for all cameras
	std::vector< Math::Vector< double, 3 > > p3dLocal;

	// local 3dpoints for each cameras (filtered - only available if observation exists)
	std::vector< std::vector < Math::Vector< double, 3 > > > p3dLocalFiltered ( numberCameras );

	// local 2d points for each camera
	std::vector< std::vector < Math::Vector< double, 2 > > > p2dLocal ( numberCameras );

	// observationCountTotal will count the number of total observations (=corners with weight != 0) in all cameras
	std::size_t observationCountTotal( 0 );
	// observations will count the number of observations (=corners with weight != 0) in each cameras
	std::vector < std::size_t > observationCount ( numberCameras );
	std::fill (observationCount.begin(), observationCount.end(), 0);


	for ( std::size_t cameraIndex = 0; cameraIndex < numberCameras; cameraIndex++ ) {
		for ( std::size_t pointIndex = startIndex; pointIndex <= endIndex; pointIndex++ ) {

			if ( points2dWeights[ cameraIndex ][ pointIndex ] != 0.0 )
			{
				OPT_LOG_TRACE( "Observation: marker corner " << pointIndex << " -> camera " << cameraIndex << ", weight=" << points2dWeights[ cameraIndex ][ pointIndex ] << ", m=" << points2d[ cameraIndex ][ pointIndex ] );
				OPT_LOG_TRACE( "According 3D Point: "<< points3d.at ( pointIndex ) );
				observations.push_back( std::make_pair( pointIndex - startIndex, cameraIndex ) );			
				p2dLocal.at(cameraIndex).push_back ( points2d.at(cameraIndex).at(pointIndex) );
				p3dLocalFiltered.at(cameraIndex).push_back ( points3d.at(pointIndex) );
				observationCountTotal++;
				observationCount.at(cameraIndex)++;
			}

			// Add 3d model points only once
			if (cameraIndex == 0) {
				p3dLocal.push_back( points3d.at( pointIndex ) ); 
			}

		}
	}

	OPT_LOG_DEBUG( observationCountTotal<<" observations found.");

	std::vector < std::size_t >::iterator minElement = std::min_element ( observationCount.begin(), observationCount.end());
	std::size_t minObs = *minElement;
	std::size_t maxObsIndex = (std::max_element (observationCount.begin(), observationCount.end())) - observationCount.begin();
	std::size_t maxObs = observationCount.at ( maxObsIndex );

	if (minObs >= minCorrespondences && (hasInitialPoseProvided || maxObs >= 4)) {
		// For the initial pose use camera with most observations
		

		// Compute initial pose
		if (!hasInitialPoseProvided) {
			OPT_LOG_DEBUG(  "Compute initial pose with "<<p2dLocal.at(maxObsIndex).size() << " observations for camera " << maxObsIndex );
			initialPose = camPoses.at( maxObsIndex ) * Algorithm::PoseEstimation2D3D::computePose( p2dLocal.at( maxObsIndex) , p3dLocalFiltered.at( maxObsIndex) ,
				camMatrices.at( maxObsIndex ), PoseEstimation2D3D::PLANAR_HOMOGRAPHY ); // there are no scoped enums in C++98 (only in C++0x onwards)
			OPT_LOG_DEBUG(  "Initial pose "<<initialPose );
		}

		// Now create the measurement vector from the local 2d points for LM optimization
		Math::Vector< double > measurements( 2 * observationCountTotal );
		std::size_t iIndex = 0;
		for ( std::size_t cameraIndex = 0; cameraIndex < numberCameras; cameraIndex++ ) {
			for ( std::size_t pointIndex = 0; pointIndex < p2dLocal.at(cameraIndex).size(); pointIndex++ ) {
				ublas::subrange( measurements, 2 * iIndex, 2 * (iIndex+1) ) = p2dLocal.at( cameraIndex ).at( pointIndex );
				OPT_LOG_TRACE( "Index: "<<iIndex << " pointIndex: "<<pointIndex);
				iIndex++;
			}
		}

		// create camera matrices and rotations for LM optimization
		std::vector< Math::Matrix< double, 3, 3 > > camRotations( camPoses.size() );
		std::vector< Math::Vector< double, 3 > > camTranslations( camPoses.size() );
		for ( std::size_t cameraIndex = 0; cameraIndex < numberCameras; cameraIndex++ )
		{
			OPT_LOG_DEBUG( "Camera "<<cameraIndex<<" pose:"  << camPoses[ cameraIndex ] );
			OPT_LOG_DEBUG( "Camera "<<cameraIndex<<" matrix:"  << camMatrices[ cameraIndex ] );

			camRotations[ cameraIndex ] = Math::Matrix< double, 3, 3 >( camPoses[ cameraIndex ].rotation() );
			camTranslations[ cameraIndex ] = camPoses[ cameraIndex ].translation();
		}

		// starting optimization
		OPT_LOG_DEBUG( "Optimizing pose over " << numberCameras << " cameras using " << observationCountTotal << " observations" );

		ObjectiveFunction< double > f( p3dLocal, camRotations, camTranslations, camMatrices, observations );
		Math::Vector< double, 6 > param;
		ublas::subrange( param, 0, 3 ) = initialPose.translation();
		ublas::subrange( param, 3, 6 ) = initialPose.rotation().toLogarithm();

		double res = Math::Optimization::levenbergMarquardt( f, param, measurements, Math::Optimization::OptTerminate( 10, 1e-6 ), Math::Optimization::OptNoNormalize() );

        // Create an error pose with covariance matrix that has the residual on its diagonal entries
		Math::ErrorPose finalPose( Math::Quaternion::fromLogarithm( ublas::subrange( param, 3, 6 ) ), ublas::subrange( param, 0, 3 ), Math::Matrix< double, 6, 6 >::identity( ) * res );
		OPT_LOG_DEBUG( "Estimated pose: " << finalPose << ", residual: " << res );

		// Everything went fine -  set weight to 1.0 and return pose
		return std::make_pair < Math::ErrorPose, double > (finalPose, res);
	} else { // not enough observations

		LOG4CPP_DEBUG( logger, "Not enough observations. Only "<<minObs<<" observations available for some camera ");
		return std::make_pair < Math::ErrorPose, double > (Math::ErrorPose(), -1.0);
	}
}

void checkConsistency (
	const std::vector < Math::Vector< double, 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector< double, 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar< double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< double, 3, 3 > >& camMatrices
	)
{
	size_t numberCameras = points2dWeights.size();

	if ( points3d.size() < 3 )
		UBITRACK_THROW( "2D6D pose estimation requires at least 3 points" );

	/*if ( points2d.size() < 2 ) {
		LOG4CPP_ERROR ( logger, "2D6D pose estimation requires at least 2 sets of 2d points. Number of sets provided: "<<points2d.size() );
		UBITRACK_THROW( "2D6D pose estimation requires at least 2 sets of 2d points" );
	}*/

	if ( points2d.size() != points2dWeights.size() || points2d.size() != camPoses.size() || points2d.size() != camMatrices.size() )
		UBITRACK_THROW( "All input sets must have the same number of cameras" );

	for ( size_t cameraIndex = 0; cameraIndex < numberCameras; cameraIndex++ )
	{
		if ( points3d.size() != points2dWeights[ cameraIndex ].size() || points3d.size() != points2d[ cameraIndex ].size() )
			UBITRACK_THROW( "All cameras must have same number of measurements as 3D points" );
	}
}


void multipleCameraPoseEstimationWithLocalBundles (
	const std::vector < Math::Vector< double, 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector< double, 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar< double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< double, 3, 3 > >& camMatrices,
	const int minCorrespondences,
	std::vector < Math::ErrorPose >& poses,
	std::vector < Math::Scalar < double > >& poseWeights,
	std::vector < Math::Scalar < int > >& localBundleSizes
	)
{
	namespace ublas = boost::numeric::ublas;
	size_t numberCameras = points2dWeights.size();

	checkConsistency ( points3d, points2d, points2dWeights, camPoses, camMatrices );

	// Offset for the current local bundle	
	int localBundleOffset = 0;

	LOG4CPP_DEBUG( logger, "Processing " << localBundleSizes.size() << " local bundles..." );

	for (int localBundleIndex = 0; localBundleIndex < localBundleSizes.size(); ++localBundleIndex) {

		// Get the number of marker corners for the current cube
		int localBundleSize = localBundleSizes.at( localBundleIndex );

		LOG4CPP_DEBUG( logger, "Local bundle " << localBundleIndex <<" has "<<localBundleSize<< " 2d points. Offset in global bundle list: "<<localBundleOffset);

		std::pair < Math::ErrorPose , double > estimate = 
			multipleCameraEstimatePose (points3d, points2d, points2dWeights, camPoses, camMatrices, minCorrespondences, false,
			Math::Pose(), localBundleOffset, localBundleOffset + localBundleSizes.at( localBundleIndex ) - 1);

		poses.push_back (estimate.first);
		poseWeights.push_back (estimate.second);

		localBundleOffset += localBundleSize;
	}
}


void multipleCameraPoseEstimation (
	const std::vector < Math::Vector< double, 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector< double, 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar< double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< double, 3, 3 > >& camMatrices,
	const int minCorrespondences,
	Math::ErrorPose& pose,
	Math::Scalar < double > & poseWeight,
	bool hasInitialPoseProvided,
	const Math::Pose initialPose
	)
{
	checkConsistency ( points3d, points2d, points2dWeights, camPoses, camMatrices );

	std::pair < Math::ErrorPose , double > estimate = 
		multipleCameraEstimatePose (points3d, points2d, points2dWeights, camPoses, camMatrices, minCorrespondences, hasInitialPoseProvided,
			initialPose);
	pose = estimate.first;
	poseWeight = estimate.second;
}

#endif // HAVE_LAPACK

} } // namespace Ubitrack::Algorithm

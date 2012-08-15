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
 * Functions for bundle adjustment.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_CALIBRATION_BUNDLEADJUSTMENT_H_INCLUDED__
#define __UBITRACK_CALIBRATION_BUNDLEADJUSTMENT_H_INCLUDED__



#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/Matrix.h>
#include <utMath/Vector.h>
#include <vector>

namespace Ubitrack { namespace Calibration {

/**
 * Describes the structure of a bundle adjustment problem.
 * The first body pose defines the world coordinate frame. Thus, there must always be at least one body!
 */
template< class T >
struct BundleAdjustmentNetwork
{
	// some typedefs to make life easier
	typedef std::vector< Math::Matrix< 3, 3, T > > IntrinsicList;
	typedef std::vector< Math::Vector< 4, T > > DistortionList;
	typedef std::vector< Math::Vector< 3, T > > PointsList;
	typedef std::vector< std::vector< Math::Vector< 3, T > > > BodyList;
	typedef std::vector< Math::Pose > PoseList;

	/**
	 * Describes a single 2D measurement of a free 3D point within the network
	 */
	struct FreePointMeasurement
	{
		FreePointMeasurement( unsigned _iPoint, unsigned _iImage, unsigned _iCamera, const Math::Vector< 2, T >& m )
			: iPoint( _iPoint )
			, iImage( _iImage )
			, iCamera( _iCamera )
			, measurement( m )
		{}

		unsigned iPoint;
		unsigned iImage;
		unsigned iCamera;
		
		Math::Vector< 2, T > measurement;
	};
	
	/**
	 * Describes a single 2D measurement of a rigid body point within the network
	 */
	struct BodyPointMeasurement
	{
		BodyPointMeasurement( unsigned _iBody, unsigned _iPoint, unsigned _iBodyPose, unsigned _iImage, 
			unsigned _iCamera, const Math::Vector< 2, T >& m )
			: iBody( _iBody )
			, iPoint( _iPoint )
			, iBodyPose( _iBodyPose )
			, iImage( _iImage )
			, iCamera( _iCamera )
			, measurement( m )
		{}

		unsigned iBody;
		unsigned iPoint;
		unsigned iBodyPose;
		unsigned iImage;
		unsigned iCamera;
		
		Math::Vector< 2, T > measurement;
	};
	
	BundleAdjustmentNetwork( PointsList& _points, PoseList& _images, IntrinsicList& _intrinsics, 
		DistortionList& _distortions, BodyList& _bodies, PoseList& _bodyPoses )
		: points( _points )
		, bodyPoses( _bodyPoses )
		, images( _images )
		, intrinsics( _intrinsics )
		, distortions( _distortions )
		, bodies( _bodies )
		, bEstimateIntrinsics( true )
	{}


	// the following members are optimized by the bundle adjustment process:
	
	/** 3D positions of the free points */
	PointsList& points;
	
	/** list of body poses. The first defines the world coordinate frame. */
	PoseList& bodyPoses;
	
 	/** list of camera poses */
	PoseList& images;
	
	/** list of camera intrinsics matrices. Note that the lower right element is assumed to be -1! */
	IntrinsicList& intrinsics;
	
	/** list of camera distortion vectors. Must have the same size as \c intrinsics */
	DistortionList& distortions;


	// the lists of measurements

	/** list of free 3d point measurements */
	std::vector< FreePointMeasurement > freePointMeasurements;

	/** list of rigid body point measurements */
	std::vector< BodyPointMeasurement > bodyPointMeasurements;


	// other parameters of the BA
	
	/** list of rigid body configurations */
	BodyList& bodies;
	
	/** should the camera intrinsics (+distortion) be estimated? */
	bool bEstimateIntrinsics;
};


/** performs the bundle adjustment */
UBITRACK_EXPORT float bundleAdjustment( BundleAdjustmentNetwork< float >& net );
UBITRACK_EXPORT double bundleAdjustment( BundleAdjustmentNetwork< double >& net );

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

#endif

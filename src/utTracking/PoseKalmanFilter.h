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
 * @ingroup tracking
 * @file
 * Class for kalman filtering of poses
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_TRACKING_POSEKALMANFILTER_H_INCLUDED__
#define __UBITRACK_TRACKING_POSEKALMANFILTER_H_INCLUDED__


#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/ErrorVector.h>
#include <utMeasurement/Measurement.h>
#include "LinearPoseMotionModel.h"

namespace Ubitrack { namespace Tracking {

/*
 * long list of todos:
 * - allow arbitrary number of derivatives
 * - include position (well, it's no longer RotationOnly then...)
 * - allow configuration of process noise (insted of hard-coded)
 * - input and output of covariances (insted of hard-coded)
 * - correct integration of measurements from the past
 */
 
/**
 * Performs fusion of pose, position, rotation, angular velocity, etc. into a single pose.
 * The motion model is based on position and rotation derivatives of arbitrary degree.
 *
 * Unfortunately, covariances of measurements and process noise can't be 
 * configured at the moment.
 */
class UBITRACK_EXPORT PoseKalmanFilter
{
public:
	/** 
	 * Constructor.
	 * @param motionModel Motion model that defines the process noise and how many derivatives for position
	 * and orientation to use.
	 * @param bInsideOut if true, a motion model is used that assumes a correlation between orientation and 
	 * translation, e.g. when a non-moveable object is tracked by a mobile camera.
	 */
	PoseKalmanFilter( const LinearPoseMotionModel& motionModel, bool bInsideOut = false );

	/** 
	 * integrate an absolute pose measurement.
	 * @param m the measured pose with timestamp and error
	 */
	void addPoseMeasurement( const Measurement::ErrorPose& m );

	/** 
	 * integrate an absolute position measurement.
	 * @param m the measured position with timestamp
	 */
	void addPositionMeasurement( const Measurement::Position& m );

	/** 
	 * integrate an absolute rotation measurement.
	 * @param m the measured rotation with timestamp
	 */
	void addRotationMeasurement( const Measurement::Rotation& m );

	/** 
	 * integrate an angular velocity measurement.
	 * @param m the measured angular velocity with timestamp
	 */
	void addRotationVelocityMeasurement( const Measurement::RotationVelocity& m );
	
	/** 
	 * Integrate an inverted angular velocity measurement.
	 * To be used for inside-out fusion.
	 * @param m the measured angular velocity with timestamp
	 */
	void addInverseRotationVelocityMeasurement( const Measurement::RotationVelocity& m );
	
	/**
	 * compute a rotation for a given time, which may lie in the future
	 */
	Measurement::ErrorPose predictPose( Measurement::Timestamp t );

	/** type of internal state representation */
	typedef Math::Vector< 0, double > StateType;

	/** type of internal state representation */
	typedef Math::Matrix< 0, 0, double > CovarianceType;

	/** returns the internal state */
	const StateType& getState() const
	{ return m_state; }

	/** returns the internal covariance */
	const CovarianceType& getCovariance() const
	{ return m_covariance; }

	/** returns the motion model */
	const LinearPoseMotionModel& getMotionModel() const
	{ return m_motionModel; }

	/** 
	 * Performs a time update of the internal state.
	 *
	 * Note: usually, there is no need to call this method, as it is implicitly called by the
	 * addXxxMeasurement methods.
	 */
	void timeUpdate( Measurement::Timestamp t );

protected:
	
	/** normalizes the state */
	void normalize();
	
	/** the motion model */
	LinearPoseMotionModel m_motionModel;

	/** inside-out motion model? */
	bool m_bInsideOut;

	/** the state */
	StateType m_state;

	/** the covariance */
	CovarianceType m_covariance;
	
	/** timestamp of the current state */
	Measurement::Timestamp m_time;
};

} } // namespace Ubitrack::Tracking

#endif // HAVE_LAPACK

#endif

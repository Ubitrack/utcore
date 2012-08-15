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
 * Class for kalman filtering of orientations
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_TRACKING_ROTATIONONLYKF_H_INCLUDED__
#define __UBITRACK_TRACKING_ROTATIONONLYKF_H_INCLUDED__


#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/ErrorVector.h>
#include <utMeasurement/Measurement.h>

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
 * Performs fusion of absolute rotation events and angular velocity events 
 * using a motion model of constant angular angular velocity.
 *
 * Can be used with absolute measurements alone, too.
 *
 * Unfortunately, covariances of measurements and process noise can't be 
 * configured at the moment.
 */
class UBITRACK_EXPORT RotationOnlyKF
{
public:
	/** constructor */
	RotationOnlyKF();

	/** 
	 * integrate an absolute rotation measurement.
	 * @param m the measured rotation with timestamp
	 */
	void addRotationMeasurement( const Measurement::Rotation& m );

	/** 
	 * integrate an angular velocity measurement.
	 * @param m the measured angular velocity with timestamp
	 */
	void addVelocityMeasurement( const Measurement::RotationVelocity& m );
	
	/**
	 * compute a rotation for a given time, which may lie in the future
	 */
	Measurement::Rotation predict( Measurement::Timestamp t );

	/** get internal state (mostly for debugging) */
	const Math::Vector< 7 >& getState() const
	{ return m_state.value; }
	
protected:
	/** performs a time update of the internal state */
	void timeUpdate( Measurement::Timestamp t );
	
	/** the state */
	Math::ErrorVector< 7 > m_state;
	
	/** timestamp of the current state */
	Measurement::Timestamp m_time;
};

} } // namespace Ubitrack::Tracking

#endif // HAVE_LAPACK

#endif

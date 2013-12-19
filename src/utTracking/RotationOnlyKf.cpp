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
 * implementation of a class for kalman filtering of orientations
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include "RotationOnlyKF.h"
#ifdef HAVE_LAPACK
 
#include <utMath/Stochastic/CovarianceTransform.h>
#include <utMath/Optimization/Function/VectorNormalize.h>
#include "Function/QuaternionTimeUpdate.h"

// get a logger
#include <log4cpp/Category.hh>
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Tracking.RotationOnlyKF" ) );
#define KALMAN_LOGGING

#include "Kalman.h"

namespace ublas = boost::numeric::ublas;

namespace Ubitrack { namespace Tracking {

RotationOnlyKF::RotationOnlyKF()
{
	// initialize state
	m_state.value = Math::Vector< double, 7 >::zeros();
	m_state.value( 3 ) = 1;
	m_state.covariance = Math::Matrix< double, 7, 7 >::identity();

	m_time = 0;
}


void RotationOnlyKF::addRotationMeasurement( const Measurement::Rotation& m )
{
	// on first update, set quaternion
	if ( m_time == 0 )
		m->toVector( m_state.value );
		
	// time update: forward filter to requested timestamp
	timeUpdate( m.time() );
	
	// create measurement as ErrorVector
	Math::ErrorVector< double, 4 > v;
	m->negateIfCloser( Math::Quaternion::fromVector( m_state.value ) ).toVector( v.value );
	v.covariance = Math::Matrix< double, 4, 4 >::identity() * 0.004; // magic number, tune here
	
	// measurement update:
	kalmanMeasurementUpdateIdentity< 7, 4 >( m_state, v, 0, 4 );

	// normalize quaternion
	Math::Stochastic::transformRangeInternalWithCovariance< 7 >( Math::Optimization::Function::VectorNormalize( 4 ), m_state, 0, 4, 0, 4 );
}


void RotationOnlyKF::addVelocityMeasurement( const Measurement::RotationVelocity& m )
{
	// time update: forward filter to requested timestamp
	timeUpdate( m.time() );
	
	// create measurement as ErrorVector
	Math::ErrorVector< double, 3 > v;
	v.value = *m;
	v.covariance = Math::Matrix< double, 3, 3 >::identity() * 0.00001; // magic number, tune here
	
	// measurement update:
	kalmanMeasurementUpdateIdentity< 7, 3 >( m_state, v, 4, 7 );

	// normalize quaternion
	Math::Stochastic::transformRangeInternalWithCovariance< 7 >( Math::Optimization::Function::VectorNormalize( 4 ), m_state, 0, 4, 0, 4 );
}


void RotationOnlyKF::timeUpdate( Measurement::Timestamp t )
{
	// only update time for the first measurement
	if ( m_time == 0 )
	{
		m_time = t;
		return;
	}
	
	// update state
	double dt = ( (long long int)( t - m_time ) ) * 1e-9;
	Math::Stochastic::transformRangeInternalWithCovariance< 7 >( Function::QuaternionTimeUpdate( dt ), m_state, 0, 4, 0, 7 );
	
	// add process noise
	// TODO: better motion model
	double fAbsNoise = 0.001 * ( dt * dt ); // tune here
	double fVelNoise = 4.0 * ( dt * dt ); // tune here
	ublas::subrange( m_state.covariance, 0, 4, 0, 4 ) += Math::Matrix< double, 4, 4 >::identity() * fAbsNoise;
	ublas::subrange( m_state.covariance, 4, 7, 4, 7 ) += Math::Matrix< double, 3, 3 >::identity() * fVelNoise;
	
	m_time = t;
}


Measurement::Rotation RotationOnlyKF::predict( Measurement::Timestamp t )
{
	// time update: forward filter to requested timestamp
	double dt = ( t - m_time ) * 1e-9;
	Math::ErrorVector< double, 4 > result = Math::Stochastic::transformWithCovariance< 4, 7 >( Function::QuaternionTimeUpdate( dt ), m_state );
	
	return Measurement::Rotation( t, Math::Quaternion::fromVector( result.value ).normalize() );
}


} } // namespace Ubitrack::Tracking

#endif // HAVE_LAPACK

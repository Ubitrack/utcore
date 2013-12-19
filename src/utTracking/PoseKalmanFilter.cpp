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
 * implementation of a class for kalman filtering of poses
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include "PoseKalmanFilter.h"
#ifdef HAVE_LAPACK
 
#include <utMath/Stochastic/CovarianceTransform.h>
#include <utMath/Optimization/Function/VectorNormalize.h>
#include <utUtil/Exception.h>
#include "Function/PoseTimeUpdate.h"
#include "Function/InsideOutPoseTimeUpdate.h"
#include "Function/InvertRotationVelocity.h"

// get a logger
#include <log4cpp/Category.hh>
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Tracking.PoseKalmanFilter" ) );

#define KALMAN_LOGGING
#include "Kalman.h"

namespace ublas = boost::numeric::ublas;

namespace {

// define the pose measurement class
struct PoseMeasurement
{
	int m_rotStart;

	PoseMeasurement( int rotStart )
		: m_rotStart( rotStart )
	{}

	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		ublas::subrange( result, 0, 3 ) = ublas::subrange( input, 0, 3 );
		ublas::subrange( result, 3, 7 ) = ublas::subrange( input, m_rotStart, m_rotStart + 4 );
		ublas::subrange( jacobian, 0, 3, 0, 3 ) = Ubitrack::Math::Matrix< double, 3, 3 >::identity();
		ublas::subrange( jacobian, 0, 3, 3, m_rotStart + 4 ) = Ubitrack::Math::Matrix< double, 0, 0 >::zeros( 3, m_rotStart + 4 - 3 );
		ublas::subrange( jacobian, 3, 7, 0, m_rotStart ) = Ubitrack::Math::Matrix< double, 0, 0 >::zeros( 4, m_rotStart );
		ublas::subrange( jacobian, 3, 7, m_rotStart, m_rotStart + 4 ) = Ubitrack::Math::Matrix< double, 4, 4 >::identity();
	}
};

}

namespace Ubitrack { namespace Tracking {

PoseKalmanFilter::PoseKalmanFilter( const LinearPoseMotionModel& motionModel, bool bInsideOut )
	: m_motionModel( motionModel )
	, m_bInsideOut( bInsideOut )
	, m_state( Math::Vector< double >::zeros( motionModel.stateSize() ) )
	, m_covariance( Math::Matrix< double, 0, 0 >::identity( motionModel.stateSize() ) )
	, m_time( 0 )
{
	if ( bInsideOut && ( m_motionModel.posOrder() > 1 || m_motionModel.oriOrder() != 1 ) )
		UBITRACK_THROW( "PoseKalmanFilter needs posOrder==1 or 0 and oriOrder==1 when inside-out mode is used!" );

	// initialize state (apart from the zeroing above)
	if ( m_motionModel.oriOrder() >= 0 )
		m_state( 6 + 3 * m_motionModel.posOrder() ) = 1;
}


void PoseKalmanFilter::addPoseMeasurement( const Measurement::ErrorPose& m )
{
	assert( m_motionModel.posOrder() >= 0 && m_motionModel.oriOrder() >= 0 );
	int iR = 3 * ( m_motionModel.posOrder() + 1 ); // shortcut for first index of orientation
	ublas::vector_range< StateType > rotSubState( m_state, ublas::range( iR, iR + 4 ) );

	// on first update, set pose
	if ( m_time == 0 )
	{
		ublas::subrange( m_state, 0, 3 ) = m->translation();
		m->rotation().toVector( rotSubState );

		// TODO: set covariance, too
	}
		
	// time update: forward filter to requested timestamp
	timeUpdate( m.time() );
	
	// create measurement as ErrorVector
	Math::ErrorVector< double, 7 > v;
	m->toAdditiveErrorVector( v );
	LOG4CPP_TRACE( logger, "Additive covariance: " << v.covariance );
	
	// negate quaternion
	if ( ublas::inner_prod( rotSubState, ublas::subrange( v.value, 3, 7 ) ) < 0 )
		ublas::subrange( v.value, 3, 7 ) *= -1;
	
	// measurement update:
	kalmanMeasurementUpdate( m_state, m_covariance, PoseMeasurement( iR ), v.value, v.covariance, 0, iR + 4 );

	// normalize quaternion
	normalize();
}


void PoseKalmanFilter::addRotationMeasurement( const Measurement::Rotation& m )
{
	assert( m_motionModel.oriOrder() >= 0 );
	int iR = 3 + 3 * m_motionModel.posOrder(); // shortcut for first index of orientation
	ublas::vector_range< StateType > rotSubState( m_state, ublas::range( iR, iR + 4 ) );
	
	// on first update, set quaternion
	if ( m_time == 0 )
		m->toVector( rotSubState );
		
	// time update: forward filter to requested timestamp
	timeUpdate( m.time() );
	
	// create measurement as ErrorVector
	Math::ErrorVector< double, 4 > v;
	m->toVector( v.value );
	v.covariance = Math::Matrix< double, 4, 4 >::identity() * 0.004; // magic number, tune here

	// invert quaternion if necessary
	if ( ublas::inner_prod( rotSubState, v.value ) < 0 )
		v.value *= -1;
	
	// measurement update:
	kalmanMeasurementUpdateIdentity( m_state, m_covariance, v.value, v.covariance, iR, iR + 4 );

	// normalize quaternion
	normalize();
}


void PoseKalmanFilter::addRotationVelocityMeasurement( const Measurement::RotationVelocity& m )
{
	assert( m_motionModel.oriOrder() >= 1 );
	
	if ( m_time == 0 )
		return;
	
	// time update: forward filter to requested timestamp
	timeUpdate( m.time() );
	
	// create measurement as ErrorVector
	Math::ErrorVector< double, 3 > v;
	v.value = *m;
	v.covariance = Math::Matrix< double, 3, 3 >::identity() * 1e-11; // magic number, tune here
	
	// measurement update:
	int iV = 4 + 3 * ( m_motionModel.posOrder() + 1 ); // shortcut for first index of rotation velocity
	kalmanMeasurementUpdateIdentity( m_state, m_covariance, v.value, v.covariance, iV, iV + 3 );

	// normalize quaternion
	normalize();
}


void PoseKalmanFilter::addInverseRotationVelocityMeasurement( const Measurement::RotationVelocity& m )
{
	assert( m_motionModel.oriOrder() >= 1 );
	
	if ( m_time == 0 )
		return;
	
	// time update: forward filter to requested timestamp
	timeUpdate( m.time() );
	
	// create measurement as ErrorVector
	Math::ErrorVector< double, 3 > v;
	v.value = *m;
	v.covariance = Math::Matrix< double, 3, 3 >::identity() * 1e-11; // magic number, tune here
	
	// measurement update:
	int iR = 3 * ( m_motionModel.posOrder() + 1 ); // shortcut for first index of orientation
	kalmanMeasurementUpdate( m_state, m_covariance, Function::InvertRotationVelocity(), v.value, v.covariance, iR, iR + 7 );

	// normalize quaternion
	normalize();
}


void PoseKalmanFilter::timeUpdate( Measurement::Timestamp t )
{
	// only update time for the first measurement
	if ( m_time == 0 || m_time == t )
	{
		m_time = t;
		return;
	}
	
	// update state
	double dt = ( (long long int)( t - m_time ) ) * 1e-9;
	LOG4CPP_DEBUG( logger, "Time update to t = " << t << ", dt = " << dt );
	if ( m_bInsideOut )
		Math::Stochastic::transformWithCovariance( 
			Function::InsideOutPoseTimeUpdate( dt, m_motionModel.posOrder() ), 
				m_state, m_covariance, m_state, m_covariance );
	else
		Math::Stochastic::transformWithCovariance( 
			Function::PoseTimeUpdate( dt, m_motionModel.posOrder(), m_motionModel.oriOrder() ),
				m_state, m_covariance, m_state, m_covariance );
	
	// add process noise
	m_motionModel.addNoise( m_covariance, dt );
	
	m_time = t;
}


void PoseKalmanFilter::normalize()
{
	int iR = 3 * ( m_motionModel.posOrder() + 1 );
	if ( m_motionModel.oriOrder() >= 0 )
	{
		// normalize quaternion
		Math::Stochastic::transformRangeInternalWithCovariance( Math::Optimization::Function::VectorNormalize( 4 ), 
			m_state, m_covariance, iR, iR + 4, iR, iR + 4 );
	}

	if ( m_motionModel.oriOrder() >= 1 )
	{
		// check if rotation velocity is too big and reset it in this case
		if ( ublas::norm_2( ublas::subrange( m_state, iR + 4, iR + 7 ) ) > 10.0 )
		{
			LOG4CPP_NOTICE( logger, "Kalman Filter orientation instability detected. Resetting orientation derivatives." );
			ublas::subrange( m_state, iR + 4, m_state.size() ) = Math::Vector< double >::zeros( 3 * m_motionModel.oriOrder() );
		}
	}
	
		// check if position values have gotten too big and reset filter then
/* 	if ( m_motionModel.posOrder() >= 0 )
		if ( fabs( m_state( 0 ) ) > 1e3 || fabs( m_state( 1 ) ) > 1e3 || fabs( m_state( 2 ) ) > 1e3 )
		{
			LOG4CPP_NOTICE( logger, "Kalman Filter position instability detected. Resetting." );
			m_state = Math::Vector< double >::zeros( m_motionModel.stateSize() );
			m_covariance = Math::Matrix< double, 0, 0 >::identity( m_motionModel.stateSize() );
			m_time = 0;
		}
 */	
}


Measurement::ErrorPose PoseKalmanFilter::predictPose( Measurement::Timestamp t )
{
	// have measurements?
	if ( !m_time )
		UBITRACK_THROW( "kalman filter not (yet) initialized" );

	assert( m_motionModel.posOrder() >= 0 && m_motionModel.oriOrder() >= 0 );
	const int iR = 3 * ( m_motionModel.posOrder() + 1 );

	double dt = ( (long long int)( t - m_time ) ) * 1e-9;
	LOG4CPP_DEBUG( logger, "predicting for t=" << t << ", dt=" << dt );

	// update state
	Math::Vector< double > newState( m_state.size() );
	Math::Matrix< double, 0, 0 > newCovariance( m_state.size(), m_state.size() );
	if ( m_bInsideOut )
		Math::Stochastic::transformWithCovariance( 
			Function::InsideOutPoseTimeUpdate( dt, m_motionModel.posOrder() ), 
				newState, newCovariance, m_state, m_covariance );
	else
		Math::Stochastic::transformWithCovariance( 
			Function::PoseTimeUpdate( dt, m_motionModel.posOrder(), m_motionModel.oriOrder() ),
				newState, newCovariance, m_state, m_covariance );
	
	// add process noise
	m_motionModel.addNoise( newCovariance, dt );
	LOG4CPP_TRACE( logger, "predicted state:" << newState );
	
	// convert to 7x7 error
	Math::ErrorVector< double, 7 > newPose;
	ublas::subrange( newPose.value, 0, 3 ) = ublas::subrange( newState, 0, 3 );
	ublas::subrange( newPose.value, 3, 7 ) = ublas::subrange( newState, iR, iR + 4 );
	ublas::subrange( newPose.covariance, 0, 3, 0, 3 ) = ublas::subrange( newCovariance, 0, 3, 0, 3 );
	ublas::subrange( newPose.covariance, 0, 3, 3, 7 ) = ublas::subrange( newCovariance, 0, 3, iR, iR + 4 );
	ublas::subrange( newPose.covariance, 3, 7, 0, 3 ) = ublas::subrange( newCovariance, iR, iR + 4, 0, 3 );
	ublas::subrange( newPose.covariance, 3, 7, 3, 7 ) = ublas::subrange( newCovariance, iR, iR + 4, iR, iR + 4 );
	
	LOG4CPP_DEBUG( logger, "predicted pose and covariance: " << newPose.value << std::endl << newPose.covariance );
	return Measurement::ErrorPose( t, Math::ErrorPose::fromAdditiveErrorVector( newPose ) );
}


} } // namespace Ubitrack::Tracking

#endif // HAVE_LAPACK

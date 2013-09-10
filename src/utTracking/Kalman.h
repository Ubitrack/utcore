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
 * Defines the measurement update of a kalman filter
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_TRACKING_KALMAN_H_INCLUDED__
#define __UBITRACK_TRACKING_KALMAN_H_INCLUDED__


#ifdef HAVE_LAPACK
 
#include <boost/numeric/bindings/lapack/posv.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <utMath/ErrorVector.h>

// to turn on logging of internal processing, create a log4cpp::Category object called "logger"
// and #define KALMAN_LOGGING before including this header 
#ifdef KALMAN_LOGGING
	#include <boost/numeric/ublas/io.hpp>
	#define KALMAN_LOG_TRACE( message ) LOG4CPP_TRACE( logger, message )
	#define KALMAN_LOG_DEBUG( message ) LOG4CPP_DEBUG( logger, message )
	#define KALMAN_LOG_NOTICE( message ) LOG4CPP_NOTICE( logger, message )
#else
	#define KALMAN_LOG_TRACE( message ) 
	#define KALMAN_LOG_DEBUG( message ) 
	#define KALMAN_LOG_NOTICE( message ) 
#endif

namespace Ubitrack { namespace Tracking {

/**
 * Heavily templated function that performs a measurement update of a kalman filter.
 *
 * @param state reference to state vector. Should contain the predicted value of a time update and 
 *     will be updated by the measurement.
 * @param stateCov reference to state vector covariance matrix.
 * @param measurementFunction function object of the measurement function. Must be modeled after 
 *   the \c Ubitrack::Calibration::Function::Prototype
 * @param measurement reference to measurement vector
 * @param measurementCov reference to measurement covariance
 * @param iBegin index of first element of state vector used as input to the measurement function (usually 0).
 * @param iEnd index of element after subvector of state used as input to measurement function (usually state.size()).
 */
template< class VState, class MStateCov, class MF, class VMeas, class MMeasCov >
void kalmanMeasurementUpdate( VState& state, MStateCov& stateCov, const MF& measurementFunction, 
	const VMeas& measurement, const MMeasCov& measurementCov, std::size_t iBegin, std::size_t iEnd )
{
	namespace ublas = boost::numeric::ublas;
	namespace blas = boost::numeric::bindings::blas;
	namespace lapack = boost::numeric::bindings::lapack;

	// some useful shortcuts
	typedef typename VState::value_type VType;
	const std::size_t inSize( iEnd - iBegin );
	const std::size_t measSize( measurement.size() );
	const std::size_t stateSize( state.size() );
	
	KALMAN_LOG_DEBUG( "state before: " << state );
	KALMAN_LOG_DEBUG( "covariance before: " << stateCov );
	
	// compute predicted measurement and jacobian
	ublas::vector< VType > predicted( measSize );
	ublas::matrix< VType, ublas::column_major > jacobian( measSize, inSize );
	measurementFunction.evaluateWithJacobian( predicted, ublas::subrange( state, iBegin, iEnd ), jacobian );
	KALMAN_LOG_DEBUG( "predicted: " << predicted );
	KALMAN_LOG_TRACE( "jacobian: " << jacobian );
	
	// compute predicted measurement error
	ublas::matrix< VType > im( ublas::prod( jacobian, ublas::subrange( stateCov, iBegin, iEnd, iBegin, iEnd ) ) );
	ublas::matrix< VType, ublas::column_major > matInv( ublas::prod( im, ublas::trans( jacobian ) ) );
	noalias( matInv ) += measurementCov;

	KALMAN_LOG_TRACE( "before inversion: " << matInv );
	if ( lapack::potrf( 'U', matInv ) == 0 )
		lapack::potri( 'U', matInv );
	else
	{
		// problem in the cholesky decomposition, try something else instead.
		KALMAN_LOG_NOTICE( "Problem in cholesky decomposition for KF. Trying something else." );
		
		// potrf has already modified matInv -> rebuild matrix
		noalias( matInv ) = ublas::prod( im, ublas::trans( jacobian ) );
		noalias( matInv ) += measurementCov;
		
		// use general matrix inversion
		ublas::vector< int > ipiv( measSize );			
		lapack::getrf( matInv, ipiv );
		lapack::getri( matInv, ipiv );
	}
	KALMAN_LOG_TRACE( "after inversion: " << matInv );
	
	// matTemp = P * H^T
	ublas::matrix< VType, ublas::column_major > matTemp( ublas::prod( 
		ublas::subrange( stateCov, 0, stateSize, iBegin, iEnd ), ublas::trans( jacobian ) ) );

	// compute kalman gain
	ublas::matrix< VType, ublas::column_major > kalmanGain( stateSize, measSize );
	blas::symm( 'R', 'U', matInv, matTemp, kalmanGain );
	KALMAN_LOG_DEBUG( "kalman gain: " << kalmanGain );

	// update state
	noalias( state ) += ublas::prod( kalmanGain, measurement - predicted );
	noalias( stateCov ) -= ublas::prod( kalmanGain, ublas::trans( matTemp ) );
	
	KALMAN_LOG_DEBUG( "state after: " << state );
	KALMAN_LOG_DEBUG( "covariance after: " << stateCov );
}


/**
 * Overload for parameters of type \c ErrorVector.
 */
template< std::size_t N, std::size_t M, class MF >
void kalmanMeasurementUpdate( Math::ErrorVector< N >& state, const MF& measurementFunction, 
	const Math::ErrorVector< M >& measurement, const std::size_t iBegin = 0, const std::size_t iEnd = N )
{
	kalmanMeasurementUpdate( state.value, state.covariance, measurementFunction, 
		measurement.value, measurement.covariance, iBegin, iEnd );
}


/**
 * measurement update for cases where the measurement function just extracts a certain sub-vector 
 * of the state, i.e. where the jacobian looks like ( 0 | I | 0 )
 *
 * @param state reference to state vector. Should contain the predicted value of a time update and 
 *     will be updated by the measurement.
 * @param stateCov reference to state vector covariance matrix.
 * @param measurement reference to measurement vector
 * @param measurementCov reference to measurement covariance
 * @param iBegin index of first element of state vector to update
 * @param iEnd index of element after subvector of state vector to update
 */
template< class VState, class MStateCov, class VMeas, class MMeasCov >
void kalmanMeasurementUpdateIdentity( VState& state, MStateCov& stateCov,  
	const VMeas& measurement, const MMeasCov& measurementCov, const std::size_t iBegin, const std::size_t iEnd )
{
	namespace ublas = boost::numeric::ublas;
	namespace blas = boost::numeric::bindings::blas;
	namespace lapack = boost::numeric::bindings::lapack;

	// some useful shortcuts
	typedef typename VState::value_type VType;
	const std::size_t measSize( measurement.size() );
	const std::size_t stateSize( state.size() );

	KALMAN_LOG_DEBUG( "state before: " << state );
	KALMAN_LOG_DEBUG( "covariance before: " << stateCov );
	
	// compute predicted measurement
	ublas::vector< VType > predicted( ublas::subrange( state, iBegin, iEnd ) );
	KALMAN_LOG_DEBUG( "predicted: " << predicted );
	
	// compute predicted measurement error
	ublas::matrix< VType, ublas::column_major > matInv( ublas::subrange( stateCov, iBegin, iEnd, iBegin, iEnd ) + measurementCov );

	KALMAN_LOG_TRACE( "before inversion: " << matInv );
	if ( lapack::potrf( 'U', matInv ) == 0 )
		lapack::potri( 'U', matInv );
	else
	{
		// problem in the cholesky decomposition, try something else instead.
		KALMAN_LOG_NOTICE( "Problem in cholesky decomposition for KF. Trying something else." );
		
		// potrf has already modified matInv -> rebuild matrix
		noalias( matInv ) = ublas::subrange( stateCov, iBegin, iEnd, iBegin, iEnd ) + measurementCov;
		
		// use general matrix inversion
		ublas::vector< int > ipiv( measSize );			
		lapack::getrf( matInv, ipiv );
		lapack::getri( matInv, ipiv );
	}
	KALMAN_LOG_TRACE( "after inversion: " << matInv );
	
	// matTemp = P * H^T
	ublas::matrix< VType, ublas::column_major > matTemp( ublas::subrange( stateCov, 0, stateSize, iBegin, iEnd ) );
	
	// compute kalman gain
	ublas::matrix< VType, ublas::column_major > kalmanGain( stateSize, measSize );
	blas::symm( 'R', 'U', matInv, matTemp, kalmanGain );
	KALMAN_LOG_DEBUG( "kalman gain: " << kalmanGain );
	
	// update state
	noalias( state ) += ublas::prod( kalmanGain, measurement - predicted );
	noalias( stateCov ) -= ublas::prod( kalmanGain, ublas::trans( matTemp ) );

	KALMAN_LOG_DEBUG( "state after: " << state );
	KALMAN_LOG_DEBUG( "covariance after: " << stateCov );
}

/**
 * Overload for parameters of type \c ErrorVector.
 */
template< std::size_t N, std::size_t M >
void kalmanMeasurementUpdateIdentity( Math::ErrorVector< N >& state, 
	const Math::ErrorVector< M >& measurement, const std::size_t iBegin, const std::size_t iEnd )
{
	kalmanMeasurementUpdateIdentity( state.value, state.covariance,  
		measurement.value, measurement.covariance, iBegin, iEnd );
}


// TODO: kalmanMeasurementUpdate without iBegin and iEnd (using optimized uBlas)

} } // namespace Ubitrack::Tracking

#undef KALMAN_LOG_DEBUG
#undef KALMAN_LOG_TRACE

#endif // HAVE_LAPACK

#endif // __UBITRACK_TRACKING_KALMAN_H_INCLUDED__


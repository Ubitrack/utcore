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
 * Implements online rotation-only hand-eye-calibration
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include <log4cpp/Category.hh>
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Calibration.OnlineRotHec" ) );
#define KALMAN_LOGGING
 
#include "OnlineRotHec.h"
#ifdef HAVE_LAPACK

#include <utMath/Optimization/Function/LinearFunction.h>
#include <utMath/Stochastic/Kalman.h>

namespace Ubitrack { namespace Calibration {

namespace ublas = boost::numeric::ublas;

static void skewMatrix( Math::Matrix< double, 3, 3 >& m, const Math::Vector< double, 3 > v )
{
	m( 0, 0 ) = 0;
	m( 0, 1 ) = -v( 2 );
	m( 0, 2 ) = v( 1 );
	m( 1, 0 ) = v( 2 );
	m( 1, 1 ) = 0;
	m( 1, 2 ) = -v( 0 );
	m( 2, 0 ) = -v( 1 );
	m( 2, 1 ) = v( 0 );
	m( 2, 2 ) = 0;
}


OnlineRotHec::OnlineRotHec()
{
	m_state.value = Math::Vector< double, 3 >( 0, 0, 0 );
	m_state.covariance = Math::Matrix< double, 3, 3 >::identity() * 1e6;
}


void OnlineRotHec::addMeasurement( const Math::Quaternion& q, const Math::Quaternion& r )
{
	// make sure the signs of both w's are equal
	const double nq = q.w() < 0 ? -1 : 1;
	const double nr = r.w() < 0 ? -1 : 1;
	
	Math::ErrorVector< double, 3 > kalmanMeasurement;
	kalmanMeasurement.value( 0 ) = r.x() * nr - q.x() * nq;
	kalmanMeasurement.value( 1 ) = r.y() * nr - q.y() * nq;
	kalmanMeasurement.value( 2 ) = r.z() * nr - q.z() * nq;
	kalmanMeasurement.covariance = Math::Matrix< double, 3, 3 >::identity();

	// do the filter update
	Math::Matrix< double, 3, 3 > h;
	skewMatrix( h, Math::Vector< double, 3 >( q.x() * nq + r.x() * nr, q.y() * nq + r.y() * nr, q.z() * nq + r.z() * nr ) );
	Math::Stochastic::kalmanMeasurementUpdate< 3, 3 >( m_state, Math::Optimization::Function::LinearFunction< 3, 3, double >( h ), 
		kalmanMeasurement, 0, m_state.value.size() );
}


Math::Quaternion OnlineRotHec::computeResult() const
{
	const double n = ublas::norm_2( m_state.value );
	const double s = 1.0 / std::sqrt( 1.0 + n * n );
	const double c = std::sqrt( 1.0 - ( n * s ) * ( n * s ) );
	return Math::Quaternion( m_state.value( 0 ) * s, m_state.value( 1 ) * s, m_state.value( 2 ) * s, c );
}

} }

#endif // HAVE_LAPACK

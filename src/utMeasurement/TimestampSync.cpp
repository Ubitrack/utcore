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
 * @ingroup datastructures
 * @file implementation of timestamp synchronization
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include "TimestampSync.h"
 
//#define DEBUG_TIMESTAMP_SYNC

#ifdef DEBUG_TIMESTAMP_SYNC
#include <iostream>
#include <iomanip>
#endif

#include <cmath>
#include <algorithm>

namespace Ubitrack { namespace Measurement {

namespace {
	// parameters to tune

	/** process noise of absolute local time (in seconds/measurement) */
	static const double g_fLocalNoise = 1e-3;

	/** process noise of gain */
	static const double g_fGainNoise = 1e-4;

	/** measurement noise (seconds/measurement) */
	static const double g_fMeasurementNoise = 3e-2;

	/** weight of measurements in outlier detection */
	static const double g_fDeviationWeight = 0.01;

} // anonymous namespace


void TimestampSync::setFrequency( double approxNativeFreq, double approxLocalFreq )
{
	m_estGain = approxLocalFreq / approxNativeFreq;
	
	// compute noise values
	m_fLocalNoise = std::pow( g_fLocalNoise * approxLocalFreq, 2 );
	m_fGainNoise = std::pow( g_fGainNoise * m_estGain, 2 );
	m_fMeasurementNoise = std::pow( g_fMeasurementNoise * approxLocalFreq, 2 );
	
	// compute initial error
	m_p1 = std::pow( approxLocalFreq, 2 );
	m_p2 = 0.0;
	m_p3 = std::pow( 1e1 * m_estGain, 2 );
	
	m_avgDeviationSquare = std::pow( 0.001 * approxLocalFreq, 2 );
}

	
Timestamp TimestampSync::convertNativeToLocal( double native, Timestamp local )
{
	// initialization on first event;
	if ( m_events == 0 )
	{
		m_lastNative = native;
		m_estLocal = local;
	}

	// time update
	double deltaN = native - m_lastNative;
	m_estLocal += static_cast< long long >( deltaN * m_estGain );
	m_p1 += 2.0 * deltaN * m_p2 + deltaN * deltaN * m_p3 + m_fLocalNoise;
	m_p2 += deltaN * m_p3;
	m_p3 += m_fGainNoise;
	
	// outlier detection
	double deltaL = double( static_cast< long long >( local - m_estLocal ) );
	double deviationSquare = deltaL * deltaL;
	double thresholdSquare = m_avgDeviationSquare * 9; // 3 sigma = 99.7% confidence (chebychev: 89% if not gaussian)
	
	if ( deviationSquare < thresholdSquare || m_outlierBudget < 0 )
	{
		// measurement update
		double k1 = m_p1 / ( m_fMeasurementNoise + m_p1 );
		double k2 = m_p2 / ( m_fMeasurementNoise + m_p1 );
		
		m_estGain += k2 * deltaL;
		m_estLocal += static_cast< long long >( k1 * deltaL );
		
		m_p1 = m_fMeasurementNoise * k1;
		m_p3 -= m_p2 * k2;
		m_p2 = m_fMeasurementNoise * k2;
		
		if ( m_events > 50 )
			m_avgDeviationSquare += ( deviationSquare - m_avgDeviationSquare ) * g_fDeviationWeight;
		
		m_outlierBudget = std::min( 40, m_outlierBudget + 1 );
	}
	else
		m_outlierBudget -= 2;

#ifdef DEBUG_TIMESTAMP_SYNC
	if ( native < 1e9  )
		std::cerr << std::setprecision( 15 ) << native << " " << local / 1000 << " " << 
			static_cast< long long >( local - m_estLocal ) / 1000 << " " << 
			m_estGain << " " << m_estLocal / 1000 << " " << 
			std::sqrt( m_avgDeviationSquare ) * 3 / 1000 << std::endl;
#endif
			
	m_lastNative = native;
	m_events++;

	return m_estLocal;
}


} } // namespace Ubitrack::Measurement

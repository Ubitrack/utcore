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
 * @file
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */


#ifndef _Ubitrack_Measurement_TimestampSyncLS_INCLUDED_
#define _Ubitrack_Measurement_TimestampSyncLS_INCLUDED_

#define DEBUG_TIMESTAMP_SYNC

#ifdef DEBUG_TIMESTAMP_SYNC
#include <iostream>
#include <iomanip>
// #include <cmath>
#endif

#include <algorithm>
#include <utCore.h>
#include <utMeasurement/Timestamp.h>

namespace Ubitrack { namespace Measurement {

// parameters to tune

/** weight of new time measurements */
static const double g_fWeight = 0.001;

/** weight of measurements in outlier detection */
static const double g_fDeviationWeight = 0.01;


/**
 * see class TimestampSync.
 * This does the same thing, but using an exponentially weighted recursive least-squares algorithm.
 */
UBITRACK_EXPORT class TimestampSyncLS
{
public:
	TimestampSyncLS()
		: m_events( 0 )
		, m_avgNative( 0.0 )
		, m_avgNativeSquare( 0.0 )
		, m_avgLocal( 0.0 )
		, m_avgLocalNative( 0.0 )
		, m_outlierBudget( 0 )
		, m_avgDeviationSquare( 1e6 * 1e6 )
	{}

	Timestamp convertNativeToLocal( double native )
	{
		return convertNativeToLocal( native, now() );
	}

	Timestamp convertNativeToLocal( double native, Timestamp local )
	{
		if ( m_events == 0 )
			m_firstLocal = local;

		double fLocal = static_cast< double >( local - m_firstLocal );
		double fExtrapolated, deviationSquare, thresholdSquare;
		double fGain, var;
		
		if ( m_events > 10 )
		{
			// extrapolate
			var =  m_avgNativeSquare - m_avgNative * m_avgNative;
			double fOffset = m_avgNativeSquare * m_avgLocal - m_avgNative * m_avgLocalNative;
			fGain = m_avgLocalNative - m_avgLocal * m_avgNative;
			fExtrapolated = ( native * fGain + fOffset ) / var;
		
			// compute threshold
			double deviation = fExtrapolated - fLocal;
			deviationSquare = deviation * deviation;
			thresholdSquare = m_avgDeviationSquare * 9;
		}
		else
		{
			fExtrapolated = fLocal;
#ifdef DEBUG_TIMESTAMP_SYNC
			fGain = 1.0; var = 1.0; // only for debug printing
#endif
		}
		
		// perform averaging
		if ( m_events < 100 || deviationSquare < thresholdSquare || m_outlierBudget < 0 )
		{
			double fWeight = m_events < 100 ? 1.0 / ( m_events + 1 ) : g_fWeight;
			m_avgNative += ( native - m_avgNative ) * fWeight;
			m_avgLocal += ( fLocal - m_avgLocal ) * fWeight;
			m_avgNativeSquare += ( native * native - m_avgNativeSquare ) * fWeight;
			m_avgLocalNative += ( fLocal * native - m_avgLocalNative ) * fWeight;

			if ( m_events > 10 )
				m_avgDeviationSquare += ( deviationSquare - m_avgDeviationSquare ) * g_fDeviationWeight;
			
			m_outlierBudget = std::min( 40, m_outlierBudget + 1 );
		}
		else
			m_outlierBudget -= 2;

#ifdef DEBUG_TIMESTAMP_SYNC
 		if ( native < 1e9  )
			std::cerr << std::setprecision( 15 ) << native << " " << local / 1000 << " " << 
				static_cast< long long >( fLocal-fExtrapolated ) / 1000 << " " << 
				fGain / var * 1e-9 << " " << ( m_firstLocal + static_cast< long long >( fExtrapolated ) ) / 1000 << " " << 
				std::sqrt( m_avgDeviationSquare ) * 3 / 1000 << std::endl;
#endif
				
		m_events++;

		// return the input value for the first 10 times
		return m_firstLocal + static_cast< Timestamp >( fExtrapolated );
	}

	unsigned getEventCount() const
	{ return m_events; }
	
protected:

	unsigned m_events;
	Timestamp m_firstLocal;
	
	// variables for recursive exponentially-weighted least squares
	double m_avgNative;
	double m_avgNativeSquare;
	double m_avgLocal;
	double m_avgLocalNative;

	// variables for outlier detection
	int m_outlierBudget;
	double m_avgDeviationSquare;
};

}}

#endif // _Ubitrack_Measurement_TimestampSyncLS_INCLUDED_

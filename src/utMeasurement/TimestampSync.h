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


#ifndef _Ubitrack_Measurement_TimestampSync_INCLUDED_
#define _Ubitrack_Measurement_TimestampSync_INCLUDED_

#include <utCore.h>
#include <utMeasurement/Timestamp.h>

namespace Ubitrack { namespace Measurement {

/**
 * Class that handles synchronization between two different clocks.
 *
 * Can be used e.g. to synchronize a sensors's internal clock to the computer's local clock. 
 * The sensors native clock is assumed to be precise whereas the local timestamp can have considerable
 * jitter when not using a real-time operating system. The shift and scaling between the clocks is
 * computed online using a kalman filter.
 */
class UBITRACK_EXPORT TimestampSync
{
public:
	/**
	 * Constructor.
	 *
	 * The timestamp synchronization needs to be initialized with approximate values of the local
	 * clock and the native clock resolution. The values need not be correct, but the order of
	 * magnitude should be.
	 *
	 * @param approxNativeFreq frequency resolution of the sensor's native clock.
	 * @param approxLocalFrequency resolution of the local clock.
	 */
	TimestampSync( double approxNativeFreq, double approxLocalFreq = 1e9 )
		: m_events( 0 )
		, m_outlierBudget( -100 ) // no outlier detection the first 100 measurements
	{
		setFrequency( approxNativeFreq, approxLocalFreq );
	}

	/** Defines the units of native and local. Values need not be precise. */
	void setFrequency( double approxNativeFreq, double approxLocalFreq = 1e9 );

	/**
	 * Add a sensor timestamp and relate it to the current system clock.
	 *
	 * @param native native sensor clock value
	 * @return the sensor time converted to a local time
	 */
	Timestamp convertNativeToLocal( double native )
	{
		return convertNativeToLocal( native, now() );
	}

	/**
	 * Add a sensor timestamp and relate it to the system clock at an arbitrary time.
	 *
	 * @param native native sensor clock value
	 * @param corresponding system clock value
	 * @return the sensor time converted to a local time
	 */
	Timestamp convertNativeToLocal( double native, Timestamp local );

	/**
	 * returns the number of timestamps processed
	 */
	unsigned getEventCount() const
	{ return m_events; }
	
protected:

	// number of treated events
	unsigned m_events;

	// time of last received native time
	double m_lastNative;
	
	// process & measurement noise (depends on timer resolution)
	double m_fLocalNoise;
	double m_fGainNoise;
	double m_fMeasurementNoise;
	
	// currently estimated local time
	Timestamp m_estLocal;
	
	// currently estimated gain
	double m_estGain;
	
	// current covariance matrix = { { p1, p2 }, { p2, p3 } }
	double m_p1;
	double m_p2;
	double m_p3;

	// variables for outlier detection
	int m_outlierBudget;
	double m_avgDeviationSquare;
};

}}

#endif // _Ubitrack_Measurement_TimestampSync_INCLUDED_

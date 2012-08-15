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

#include "Timestamp.h"

#ifdef _MSC_VER

	#include <sys/timeb.h>
	#include <time.h>

	#include <utUtil/OS.h>
	#include "TimestampSync.h"

	#include <log4cpp/Category.hh>
	static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Measurement.Timestamp" ) );

#else
    #include <stdio.h>

	#include <sys/time.h>

#endif

namespace Ubitrack { namespace Measurement {


Timestamp now()
{
#ifdef _MSC_VER

	// Windows time has bad resolution (10ms), therefore we prefer to use the High Performance Counter 
	// and synchronize it to the real time clock
	
	// read HPC
	long long hiPerf = Util::getHighPerformanceCounter();
	
	// read RTC
	struct _timeb wintime;
	#if _MSC_VER < 1400
		_ftime( &wintime );
	#else
		_ftime_s( &wintime );
	#endif
	Timestamp rtc = Timestamp( wintime.time ) * 1000000000 + Timestamp( wintime.millitm ) * 1000000;

	// synchronization
	static TimestampSync synchronizer( 1e9 );
	static bool bUseHpc = true;
	static double lastHpcFreq;
	
	if ( bUseHpc )
	{
		// check if HPC is any good (has constant frequency)
		double hpcFreq = Util::getHighPerformanceFrequency();
		if ( synchronizer.getEventCount() > 0 )
			if ( hpcFreq != lastHpcFreq )
			{
				bUseHpc = false;
				LOG4CPP_WARN( logger, "Your CPU frequency is not constant (power save mode?). Timestamps will be unprecise." );
				return rtc;
			}
			else
			{}
		else
		{
			// initialization
			synchronizer.setFrequency( hpcFreq );
			lastHpcFreq = hpcFreq;
		}
		
		// convert hpc to rtc time
		return synchronizer.convertNativeToLocal( double( hiPerf ), rtc );
		//return Timestamp( hiPerf / hpcFreq * 1e9 );
	}
	else
		return rtc;
	
#else

	// Unix time
	struct timeval tv;
	gettimeofday( &tv, 0 );
	return ((Timestamp)(tv.tv_sec)) * 1000000000 + ((Timestamp)(tv.tv_usec)) * 1000;

#endif	
}


std::string timestampToString( Timestamp t )
{
	std::string result;
	time_t mytime = (time_t)(t/1000000000);
	#if !defined( _MSC_VER ) || _MSC_VER >= 1400
		char buffer[50]; // asctime requires at least 26 bytes
		struct tm temp;
	#endif

	#ifdef _MSC_VER
		#if _MSC_VER >= 1400
			gmtime_s( &temp, &mytime );
			asctime_s( buffer, sizeof(buffer), &temp );
			result = buffer;
		#else
			result = asctime( gmtime( &mytime ) );
		#endif			
	#else
		gmtime_r( &mytime, &temp );
		asctime_r( &temp, buffer );
		result = buffer;
	#endif

	result.replace( result.length()-1, 4, " UTC" );
	return result;
}


std::string timestampToShortString( Timestamp t )
{
	time_t mytime = (time_t)(t/1000000000);
	struct tm temp;
	struct tm* pTm = &temp;

	#ifdef _MSC_VER
		#if _MSC_VER >= 1400
			localtime_s( &temp, &mytime );
		#else
			pTm = localtime( &mytime );
		#endif			
	#else
		localtime_r( &mytime, &temp );
	#endif

	char buffer[ 32 ];
	sprintf( buffer, "%02d:%02d:%02d.%06d", pTm->tm_hour, pTm->tm_min, pTm->tm_sec, int( ( t / 1000 ) % 1000000 ) );
	return std::string( buffer );
}

} } // namespace Ubitrack::Measurement



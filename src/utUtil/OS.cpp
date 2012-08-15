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
 * @file
 * Implements abstractions for operating system specific functions
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#include "OS.h"

#ifdef WIN32
#include "CleanWindows.h"
#else
#include <time.h>
#include <sys/time.h>
#endif

namespace Ubitrack { namespace Util {

#ifdef _WIN32

void sleep( unsigned ms, unsigned )
{
	Sleep( ms );
}


long long getHighPerformanceCounter()
{
	LARGE_INTEGER counter;
	QueryPerformanceCounter( &counter );
	return counter.QuadPart;
}


double getHighPerformanceFrequency()
{
	LARGE_INTEGER counter;
	QueryPerformanceFrequency( &counter );
	return static_cast< double >( counter.QuadPart );
}

#else // unix

void sleep( unsigned ms, unsigned ns )
{	
	timespec t;
	t.tv_sec = ms / 1000;
	t.tv_nsec = ( ms % 1000 ) * 1000000 + ns;
	nanosleep( &t, NULL );
}


long long getHighPerformanceCounter()
{
	// TODO: read CPU counter
	struct timeval tv;
	gettimeofday( &tv, 0 );
	return tv.tv_sec * static_cast< long long >( 1000000 ) + tv.tv_usec;
}


double getHighPerformanceFrequency()
{
	return 1000000.0;
}
	
#endif

} } // namespace Ubitrack::Util

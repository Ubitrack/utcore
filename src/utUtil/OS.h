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
 * Abstractions for operating system specific functions
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 
 
#ifndef __UBITRACK_UTIL_OS_H_INCLUDED__
#define __UBITRACK_UTIL_OS_H_INCLUDED__

#include <utCore.h>

namespace Ubitrack { namespace Util {

/** 
 * Holds execution for the given number of milliseconds.
 * boost::thread::sleep sucks, that's why we are providing our own version.
 *
 * @param ms number of milliseconds to sleep
 * @param ns number of additional nanoseconds to sleep (in case finer-grained sleeping should be necessary)
 */
UBITRACK_EXPORT void sleep( unsigned ms, unsigned ns = 0 );

/** retrieves the high performance counter value */
UBITRACK_EXPORT long long getHighPerformanceCounter();

/** retrieves the high performance counter frequency in Hz */
UBITRACK_EXPORT double getHighPerformanceFrequency();

} } // namespace Ubitrack::Util

#endif

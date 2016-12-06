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
 * @ingroup dataflow_framework
 * @file
 * Windows Event Tracing Support
 *
 * @author Ulrich Eck <ueck@net-labs.de>
 * based on work from: https://randomascii.wordpress.com/2015/04/14/uiforetw-windows-performance-made-easier/
 */ 
#ifndef UBITRACK_ETWTRACING_H
#define UBITRACK_ETWTRACING_H


#include <utCore.h>

// ETW (Event Tracing for Windows) profiling helpers.
// This allows easy insertion of Generic Event markers into ETW/xperf tracing
// which then aids in analyzing the traces and finding performance problems.
// The usage patterns are to use ETWBegin and ETWEnd (typically through the
// convenience class CETWScope) to bracket time-consuming operations. In addition
// ETWFrameMark marks the beginning of each frame, and ETWMark can be used to
// mark other notable events. More event types and providers can be added as needed.


#if defined( _MSC_VER )
#pragma once
#endif

typedef long long int64;

#ifdef	_WIN32

#include <SDKDDKVer.h>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

// ETW support should be compiled in for all Windows PC platforms. It isn't
// supported on Windows XP but that is determined at run-time. This #define
// is used to let the code compile (but do nothing) on other operating systems.
#ifdef ENABLE_EVENT_TRACING
#if (WINVER >= _WIN32_WINNT_VISTA)
#define	ETW_MARKS_ENABLED
#endif
#endif

#endif  // _WIN32

#ifdef	ETW_MARKS_ENABLED

#include <sal.h> // For _Printf_format_string_

UBITRACK_EXPORT int64 __cdecl ETWUbitrackEventQueueDispatchBegin(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName);
UBITRACK_EXPORT int64 __cdecl ETWUbitrackEventQueueDispatchEnd(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName, int64 startTime);
UBITRACK_EXPORT void __cdecl ETWUbitrackEventQueueDispatchDiscard(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName);
//UBITRACK_EXPORT void __cdecl ETWUbitrackEventQueueApplication(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName,  _In_z_ PCSTR text);

UBITRACK_EXPORT void __cdecl ETWUbitrackMeasurementCreate(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName, _In_z_ PCSTR portName);
UBITRACK_EXPORT void __cdecl ETWUbitrackMeasurementReceive(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName, _In_z_ PCSTR portName);
UBITRACK_EXPORT void __cdecl ETWUbitrackAllocateCpu(unsigned int bytes);
UBITRACK_EXPORT void __cdecl ETWUbitrackAllocateGpu(unsigned int bytes);
UBITRACK_EXPORT void __cdecl ETWUbitrackGpuUpload(unsigned int bytes);
UBITRACK_EXPORT void __cdecl ETWUbitrackGpuDownload(unsigned int bytes);

#else // ETW_MARKS_ENABLED

UBITRACK_EXPORT int64 __cdecl ETWUbitrackEventQueueDispatchBegin(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName) { return 0; };
UBITRACK_EXPORT int64 __cdecl ETWUbitrackEventQueueDispatchEnd(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName, int64 startTime) { return 0; };
UBITRACK_EXPORT void __cdecl ETWUbitrackEventQueueDispatchDiscard(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName) {};
//UBITRACK_EXPORT void __cdecl ETWUbitrackEventQueueApplication(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName,  _In_z_ PCSTR text) {};

UBITRACK_EXPORT void __cdecl ETWUbitrackMeasurementCreate(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName, _In_z_ PCSTR portName) {};
UBITRACK_EXPORT void __cdecl ETWUbitrackMeasurementReceive(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName, _In_z_ PCSTR portName) {};
UBITRACK_EXPORT void __cdecl ETWUbitrackAllocateCpu(unsigned int bytes) {};
UBITRACK_EXPORT void __cdecl ETWUbitrackAllocateGpu(unsigned int bytes) {};
UBITRACK_EXPORT void __cdecl ETWUbitrackGpuUpload(unsigned int bytes) {};
UBITRACK_EXPORT void __cdecl ETWUbitrackGpuDownload(unsigned int bytes) {};


#endif // ETW_MARKS_ENABLED


#endif // UBITRACK_ETWTRACING_H
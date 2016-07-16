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


#ifdef	_WIN32
#include <stdio.h>
#include "ETWTracing.h"
#endif

#ifdef ETW_MARKS_ENABLED

#include <SDKDDKVer.h>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

// These are defined in evntrace.h but you need a Vista+ Windows
// SDK to have them available, so I define them here.
#define EVENT_CONTROL_CODE_DISABLE_PROVIDER 0
#define EVENT_CONTROL_CODE_ENABLE_PROVIDER  1
#define EVENT_CONTROL_CODE_CAPTURE_STATE    2

#ifdef EVNTAPI
// Make sure we don't include evntprov.h before #defining EVNTAPI
// If we do then we may end up calling the imported NTDLL functions and
// then our code will fail to run on Windows XP.
#error "We are calling the imported functions. This will not work on Windows XP"
#endif

// EVNTAPI is used in evntprov.h which is included by ETWProviderGenerated.h
// We define EVNTAPI without the DECLSPEC_IMPORT specifier so that
// we can implement these functions locally instead of using the import library,
// and can therefore still run on Windows XP.
#define EVNTAPI __stdcall
// Include the event register/write/unregister macros compiled from the manifest file.
// Note that this includes evntprov.h which requires a Vista+ Windows SDK.
#include "utUtil/probes_ubitrack_etw.h"

// Typedefs for use with GetProcAddress
typedef ULONG (__stdcall *tEventRegister)( LPCGUID ProviderId, PENABLECALLBACK EnableCallback, PVOID CallbackContext, PREGHANDLE RegHandle);
typedef ULONG (__stdcall *tEventWrite)( REGHANDLE RegHandle, PCEVENT_DESCRIPTOR EventDescriptor, ULONG UserDataCount, PEVENT_DATA_DESCRIPTOR UserData);
typedef ULONG (__stdcall *tEventUnregister)( REGHANDLE RegHandle );


// Helper class to dynamically load Advapi32.dll, find the ETW functions, 
// register the providers if possible, and get the performance counter frequency.
class CETWRegister
{
public:
	CETWRegister()
	{
		QueryPerformanceFrequency( &m_frequency );

		// Find Advapi32.dll. This should always succeed.
		HMODULE pAdvapiDLL = LoadLibraryW( L"Advapi32.dll" );
		if ( pAdvapiDLL )
		{
			// Try to find the ETW functions. This will fail on XP.
			m_pEventRegister = ( tEventRegister )GetProcAddress( pAdvapiDLL, "EventRegister" );
			m_pEventWrite = ( tEventWrite )GetProcAddress( pAdvapiDLL, "EventWrite" );
			m_pEventUnregister = ( tEventUnregister )GetProcAddress( pAdvapiDLL, "EventUnregister" );

			// Register our ETW providers. If registration fails then the event logging calls will fail.
			// On XP these calls will do nothing.
			// On Vista and above, if these providers have been enabled by xperf or logman then
			// the *Context globals will be modified
			// like this:
			//     MatchAnyKeyword: 0xffffffffffffffff
			//     IsEnabled: 1
			//     Level: 255
			// In other words, fully enabled.

			EventRegisterUbitrack();

			// Emit the thread ID for the main thread. This also indicates that
			// the main provider is initialized.
			//EventWriteThread_ID( GetCurrentThreadId(), "Main thread" );
		}
	}
	~CETWRegister()
	{
		// Unregister our providers.
		EventUnregisterUbitrack();
	}

	tEventRegister m_pEventRegister;
	tEventWrite m_pEventWrite;
	tEventUnregister m_pEventUnregister;

	// QPC frequency
	LARGE_INTEGER m_frequency;

} g_ETWRegister;



// Redirector function for EventRegister. Called by macros in ETWProviderGenerated.h
ULONG EVNTAPI EventRegister( LPCGUID ProviderId, PENABLECALLBACK EnableCallback, PVOID CallbackContext, PREGHANDLE RegHandle )
{
	if ( g_ETWRegister.m_pEventRegister )
		return g_ETWRegister.m_pEventRegister( ProviderId, EnableCallback, CallbackContext, RegHandle );

	return 0;
}

// Redirector function for EventWrite. Called by macros in ETWProviderGenerated.h
ULONG EVNTAPI EventWrite( REGHANDLE RegHandle, PCEVENT_DESCRIPTOR EventDescriptor, ULONG UserDataCount, PEVENT_DATA_DESCRIPTOR UserData )
{
	if ( g_ETWRegister.m_pEventWrite )
		return g_ETWRegister.m_pEventWrite( RegHandle, EventDescriptor, UserDataCount, UserData );
	return 0;
}

// Redirector function for EventUnregister. Called by macros in ETWProviderGenerated.h
ULONG EVNTAPI EventUnregister( REGHANDLE RegHandle )
{
	if ( g_ETWRegister.m_pEventUnregister )
		return g_ETWRegister.m_pEventUnregister( RegHandle );
	return 0;
}

// Call QueryPerformanceCounter
static int64 GetQPCTime()
{
	LARGE_INTEGER time;

	QueryPerformanceCounter( &time );
	return time.QuadPart;
}

// Convert a QueryPerformanceCounter delta into milliseconds
static float QPCToMS( int64 nDelta )
{
	// Convert from a QPC delta to seconds.
	float flSeconds = ( float )( nDelta / double( g_ETWRegister.m_frequency.QuadPart ) );

	// Convert from seconds to milliseconds
	return flSeconds * 1000;
}

int64 ETWUbitrackEventQueueDispatchBegin(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName)
{
	// If we are running on Windows XP or if our providers have not been enabled
	// (by xperf or other) then this will be false and we can early out.
	// Be sure to check the appropriate context for the event. This is only
	// worth checking if there is some cost beyond the EventWrite that we can
	// avoid -- the redirectors in this file guarantee that EventWrite is always
	// safe to call.
	// In this case we also avoid the potentially unreliable TLS implementation
	// (for dynamically loaded DLLs) on Windows XP.
	if ( !UBITRACK_Context.IsEnabled )
	{
		return 0;
	}

	int64 nTime = GetQPCTime();
	EventWriteEventQueueDispatchBegin( eventDomain, priority, componentName, portName );
	return nTime;
}

int64 ETWUbitrackEventQueueDispatchEnd(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName, int64 nStartTime)
{
	// If we are running on Windows XP or if our providers have not been enabled
	// (by xperf or other) then this will be false and we can early out.
	// Be sure to check the appropriate context for the event. This is only
	// worth checking if there is some cost beyond the EventWrite that we can
	// avoid -- the redirectors in this file guarantee that EventWrite is always
	// safe to call.
	// In this case we also avoid the potentially unreliable TLS implementation
	// (for dynamically loaded DLLs) on Windows XP.
	if ( !UBITRACK_Context.IsEnabled )
	{
		return 0;
	}

	int64 nTime = GetQPCTime();
	EventWriteEventQueueDispatchEnd( eventDomain, priority, componentName, portName, QPCToMS( nTime - nStartTime ) );
	return nTime;
}

void ETWUbitrackEventQueueDispatchDiscard(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName)
{
	// If we are running on Windows XP or if our providers have not been enabled
	// (by xperf or other) then this will be false and we can early out.
	// Be sure to check the appropriate context for the event. This is only
	// worth checking if there is some cost beyond the EventWrite that we can
	// avoid -- the redirectors in this file guarantee that EventWrite is always
	// safe to call.
	// In this case we also avoid the potentially unreliable TLS implementation
	// (for dynamically loaded DLLs) on Windows XP.
	if ( !UBITRACK_Context.IsEnabled )
	{
		return;
	}

	EventWriteEventQueueDispatchDiscard( eventDomain, priority, componentName, portName );
}

//void ETWUbitrackEventQueueApplication(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName,  _In_z_ PCSTR portName,  _In_z_ PCSTR text)
//{
//	// If we are running on Windows XP or if our providers have not been enabled
//	// (by xperf or other) then this will be false and we can early out.
//	// Be sure to check the appropriate context for the event. This is only
//	// worth checking if there is some cost beyond the EventWrite that we can
//	// avoid -- the redirectors in this file guarantee that EventWrite is always
//	// safe to call.
//	// In this case we also avoid the potentially unreliable TLS implementation
//	// (for dynamically loaded DLLs) on Windows XP.
//	if ( !UBITRACK_Context.IsEnabled )
//	{
//		return;
//	}
//
//	EventWriteEventQueueApplication( eventDomain, priority, componentName, portName, text );
//}


void ETWUbitrackMeasurementCreate(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName, _In_z_ PCSTR portName) {
	if (!UBITRACK_Context.IsEnabled)
	{
		return;
	}

	ETWUbitrackMeasurementCreate(eventDomain, priority, componentName, portName);
};

void ETWUbitrackMeasurementReceive(int eventDomain, unsigned long long int priority, _In_z_ PCSTR componentName, _In_z_ PCSTR portName) {
	if (!UBITRACK_Context.IsEnabled)
	{
		return;
	}

	ETWUbitrackMeasurementReceive(eventDomain, priority, componentName, portName);
};

void ETWUbitrackAllocateCpu(unsigned int bytes) {
	if (!UBITRACK_Context.IsEnabled)
	{
		return;
	}

	ETWUbitrackAllocateCpu(bytes);
};

void ETWUbitrackAllocateGpu(unsigned int bytes) {
	if (!UBITRACK_Context.IsEnabled)
	{
		return;
	}

	ETWUbitrackAllocateGpu(bytes);
};

void ETWUbitrackGpuUpload(unsigned int bytes) {
	if (!UBITRACK_Context.IsEnabled)
	{
		return;
	}

	ETWUbitrackGpuUpload(bytes);
};

void ETWUbitrackGpuDownload(unsigned int bytes) {
	if (!UBITRACK_Context.IsEnabled)
	{
		return;
	}

	ETWUbitrackGpuDownload(bytes);
};

#endif // ETW_MARKS_ENABLED
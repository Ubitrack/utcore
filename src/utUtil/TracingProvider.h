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
 * The main include for tracing the eventqueue activity
 *
 * @author Ulrich Eck <ueck@net-labs.de>
 */ 

#ifndef UBITACK_TRACINGPROVIDER_H
#define UBITACK_TRACINGPROVIDER_H

#ifdef ENABLE_EVENT_TRACING

#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#include "utUtil/probes_ubitrack_dtrace.h"
#else
//#include "utUtil/probes_ubitrack_nodtrace.h"
#endif

#ifdef HAVE_ETW
#include "ETWTracing.h"
#else
// xxx
#endif




//
//
//
//#ifdef ENABLE_EVENT_TRACING
//#ifdef HAVE_DTRACE
//			if (UBITRACK_EVENTQUEUE_DISPATCH_BEGIN_ENABLED() && pReceiverInfo ) {
//				UBITRACK_EVENTQUEUE_DISPATCH_BEGIN(m_eventDomain,
//												   messagePriority,
//												   pReceiverInfo->pPort->getComponent().getName().c_str(),
//												   pReceiverInfo->pPort->getName().c_str());
//			}
//#endif
//
//#ifdef HAVE_ETW
//			int64 _startTime = ETWUbitrackEventQueueDispatchBegin(m_eventDomain, messagePriority,
//																  pReceiverInfo->pPort->getComponent().getName().c_str(),
//																  pReceiverInfo->pPort->getName().c_str());
//#endif
//#endif
//#ifdef ENABLE_EVENT_TRACING
//#ifdef HAVE_DTRACE
//			if (UBITRACK_EVENTQUEUE_DISPATCH_BEGIN_ENABLED() && pReceiverInfo ) {
//				UBITRACK_EVENTQUEUE_DISPATCH_BEGIN(m_eventDomain,
//												   messagePriority,
//												   pReceiverInfo->pPort->getComponent().getName().c_str(),
//												   pReceiverInfo->pPort->getName().c_str());
//			}
//#endif
//
//#ifdef HAVE_ETW
//			int64 _startTime = ETWUbitrackEventQueueDispatchBegin(m_eventDomain, messagePriority,
//																  pReceiverInfo->pPort->getComponent().getName().c_str(),
//																  pReceiverInfo->pPort->getName().c_str());
//#endif
//#endif
//
//
//
//
//
//#ifdef ENABLE_EVENT_TRACING
//#ifdef HAVE_DTRACE
//			if (UBITRACK_EVENTQUEUE_DISPATCH_END_ENABLED() && pReceiverInfo ) {
//				UBITRACK_EVENTQUEUE_DISPATCH_END(m_eventDomain,
//												 messagePriority,
//												 pReceiverInfo->pPort->getComponent().getName().c_str(),
//												 pReceiverInfo->pPort->getName().c_str());
//			}
//#endif
//
//#ifdef HAVE_ETW
//			ETWUbitrackEventQueueDispatchEnd(m_eventDomain, messagePriority,
//											 pReceiverInfo->pPort->getComponent().getName().c_str(),
//											 pReceiverInfo->pPort->getName().c_str(),
//											 _startTime);
//#endif
//#endif
//#ifdef ENABLE_EVENT_TRACING
//#ifdef HAVE_DTRACE
//			if (UBITRACK_EVENTQUEUE_DISPATCH_END_ENABLED() && pReceiverInfo ) {
//				UBITRACK_EVENTQUEUE_DISPATCH_END(m_eventDomain,
//												 messagePriority,
//												 pReceiverInfo->pPort->getComponent().getName().c_str(),
//												 pReceiverInfo->pPort->getName().c_str());
//			}
//#endif
//
//
//#ifdef HAVE_ETW
//			ETWUbitrackEventQueueDispatchEnd(m_eventDomain, messagePriority,
//											 pReceiverInfo->pPort->getComponent().getName().c_str(),
//											 pReceiverInfo->pPort->getName().c_str(),
//											 _startTime);
//#endif
//#endif
//
//
//
//#ifdef ENABLE_EVENT_TRACING
//		#ifdef HAVE_DTRACE
//			if (UBITRACK_MEASUREMENT_CREATE_ENABLED()) {
//				UBITRACK_MEASUREMENT_CREATE(getEventDomain(),
//											evt.time(),
//											getName().c_str(),
//											"Output");
//			}
//#endif
//
//#ifdef HAVE_ETW
//			ETWUbitrackMeasurementCreate(getEventDomain(), evt.time(),
//										 getName().c_str(),
//										 "Output");
//#endif
//#endif
//
//
//#ifdef ENABLE_EVENT_TRACING
//			#ifdef HAVE_DTRACE
//			if (UBITRACK_MEASUREMENT_RECEIVE_ENABLED()) {
//				UBITRACK_MEASUREMENT_RECEIVE(getEventDomain(),
//											 m.time(),
//											 getName().c_str(),
//											 "Input");
//			}
//#endif
//
//#ifdef HAVE_ETW
//			ETWUbitrackMeasurementReceive(getEventDomain(), evt.time(),
//										  getName().c_str(),
//										  "Input");
//#endif
//#endif
//



#endif // ENABLE_EVENT_TRACING

#endif //UBITACK_TRACINGPROVIDER_H

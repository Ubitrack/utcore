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

#ifndef UBITRACK_TRACINGPROVIDER_H
#define UBITRACK_TRACINGPROVIDER_H

#ifdef ENABLE_EVENT_TRACING

#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#include "utUtil/probes_ubitrack_dtrace.h"
#endif


#ifdef HAVE_LTTNGUST
#include "LTTNGTracingProvider.h"
#endif

#ifdef HAVE_ETW
#include "ETWTracing.h"
#endif

#ifdef HAVE_USDT
// copied from folly/tracingin/StaticTracepoint.h
#if defined(__ELF__) &&                                                        \
    (defined(__powerpc64__) || defined(__powerpc__) || defined(__aarch64__) || \
     defined(__x86_64__) || defined(__i386__))
#include <utUtil/StaticTracepoint-ELF.h>

#define FOLLY_SDT(provider, name, ...)                                         \
  FOLLY_SDT_PROBE_N(                                                           \
    provider, name, FOLLY_SDT_NARG(0, ##__VA_ARGS__), ##__VA_ARGS__)
#else
#define FOLLY_SDT(provider, name, ...) do {} while(0)
#endif
#endif // HAVE_USDT




/*
 * TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_BEGIN(event_domain, event_priority, component_name, component_port)
 * traces the begin of the execution of an eventqueue item
 * parameters:
 * - event_domain: Which domain the item belongs to
 * - event_priority: The timestamp of the event or some other globally known identifier
 * - component_name: The component that causes the execution of the event (the receiving component of a push message)
 * - component_port: The port on which the message was received
 *
 * example:
 * TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_BEGIN(m_eventDomain, messagePriority,
 *   pReceiverInfo->pPort->getComponent().getName().c_str(), pReceiverInfo->pPort->getName().c_str())
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_BEGIN(event_domain, event_priority, component_name, component_port)\
  if (UBITRACK_EVENTQUEUE_DISPATCH_BEGIN_ENABLED()) {\
    UBITRACK_EVENTQUEUE_DISPATCH_BEGIN(event_domain, event_priority, component_name, component_port);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_BEGIN(event_domain, event_priority, component_name, component_port)\
  ___ubitrack_tracing_startTime = ETWUbitrackEventQueueDispatchBegin(event_domain, event_priority, component_name, component_port);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_BEGIN(event_domain, event_priority, component_name, component_port)\
  tracepoint(ubitrack, eventqueue_dispatch_begin, event_domain, event_priority, component_name, component_port);
#endif

#ifdef HAVE_USDT
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_BEGIN(event_domain, event_priority, component_name, component_port)\
  FOLLY_SDT(ubitrack, eventqueue_dispatch_begin, event_domain, event_priority, component_name, component_port);
#endif


/*
 * TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_END(event_domain, event_priority, component_name, component_port)
 * traces the end of the execution of an eventqueue item
 * parameters:
 * - event_domain: Which domain the item belongs to
 * - event_priority: The timestamp of the event or some other globally known identifier
 * - component_name: The component that causes the execution of the event (the receiving component of a push message)
 * - component_port: The port on which the message was received
 *
 * example:
 * TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_END(m_eventDomain, messagePriority,
 *   pReceiverInfo->pPort->getComponent().getName().c_str(), pReceiverInfo->pPort->getName().c_str())
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_END(event_domain, event_priority, component_name, component_port)\
  if (UBITRACK_EVENTQUEUE_DISPATCH_END_ENABLED()) {\
    UBITRACK_EVENTQUEUE_DISPATCH_END(event_domain, event_priority, component_name, component_port);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_END(event_domain, event_priority, component_name, component_port)\
  ETWUbitrackEventQueueDispatchEnd(event_domain, event_priority, component_name, component_port,___ubitrack_tracing_startTime);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_END(event_domain, event_priority, component_name, component_port)\
  tracepoint(ubitrack, eventqueue_dispatch_end, event_domain, event_priority, component_name, component_port);
#endif

#ifdef HAVE_USDT
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_END(event_domain, event_priority, component_name, component_port)\
  FOLLY_SDT(ubitrack, eventqueue_dispatch_end, event_domain, event_priority, component_name, component_port);
#endif

/*
 * TRACEPOINT_MEASUREMENT_CREATE(event_domain, event_priority, component_name, component_port)
 * traces the creation of a measurement
 * parameters:
 * - event_domain: Which domain the item belongs to
 * - event_priority: The timestamp of the event or some other globally known identifier
 * - component_name: The component that created the event (the sending component of a push message)
 * - component_port: The port on which the message will be sent
 *
 * example:
 * TRACEPOINT_MEASUREMENT_CREATE(getEventDomain(), evt.time(), getName().c_str(), "Output")
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_MEASUREMENT_CREATE(event_domain, event_priority, component_name, component_port)\
  if (UBITRACK_MEASUREMENT_CREATE_ENABLED()) {\
    UBITRACK_MEASUREMENT_CREATE(event_domain, event_priority, component_name, component_port);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_MEASUREMENT_CREATE(event_domain, event_priority, component_name, component_port)\
  ETWUbitrackMeasurementCreate(event_domain, event_priority, component_name, component_port);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_MEASUREMENT_CREATE(event_domain, event_priority, component_name, component_port)\
  tracepoint(ubitrack, measurement_create, event_domain, event_priority, component_name, component_port);
#endif

#ifdef HAVE_USDT
#define TRACEPOINT_MEASUREMENT_CREATE(event_domain, event_priority, component_name, component_port)\
  FOLLY_SDT(ubitrack, measurement_create, event_domain, event_priority, component_name, component_port);
#endif

/*
 * TRACEPOINT_MEASUREMENT_RECEIVE(event_domain, event_priority, component_name, component_port)
 * traces the reception of a measurement
 * parameters:
 * - event_domain: Which domain the item belongs to
 * - event_priority: The timestamp of the event or some other globally known identifier
 * - component_name: The component that receives the event (typically an application or network sink)
 * - component_port: The port on which the message was received
 *
 * example:
 * TRACEPOINT_MEASUREMENT_RECEIVE(getEventDomain(), evt.time(), getName().c_str(), "Output")
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_MEASUREMENT_RECEIVE(event_domain, event_priority, component_name, component_port)\
  if (UBITRACK_MEASUREMENT_RECEIVE_ENABLED()) {\
    UBITRACK_MEASUREMENT_RECEIVE(event_domain, event_priority, component_name, component_port);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_MEASUREMENT_RECEIVE(event_domain, event_priority, component_name, component_port)\
  ETWUbitrackMeasurementReceive(event_domain, event_priority, component_name, component_port);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_MEASUREMENT_RECEIVE(event_domain, event_priority, component_name, component_port)\
  tracepoint(ubitrack, measurement_receive, event_domain, event_priority, component_name, component_port);
#endif

#ifdef HAVE_USDT
#define TRACEPOINT_MEASUREMENT_RECEIVE(event_domain, event_priority, component_name, component_port)\
  FOLLY_SDT(ubitrack, measurement_receive, event_domain, event_priority, component_name, component_port);
#endif

 /*
 * TRACEPOINT_VISION_ALLOCATE_CPU(bytes)
 * traces the allocation on CPU memory of an image
 * parameters:
 * - bytes: number of bytes that were allocated
 *
 * example:
 * TRACEPOINT_VISION_ALLOCATE_CPU(m_width*m_height*m_channels)
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_VISION_ALLOCATE_CPU(bytes)\
  if (UBITRACK_VISION_ALLOCATE_CPU_ENABLED()) {\
    UBITRACK_VISION_ALLOCATE_CPU(bytes);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_VISION_ALLOCATE_CPU(bytes)\
  ETWUbitrackAllocateCpu(bytes);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_VISION_ALLOCATE_CPU(bytes)\
  tracepoint(ubitrack, vision_allocate_cpu, bytes);
#endif

#ifdef HAVE_USDT
#define TRACEPOINT_VISION_ALLOCATE_CPU(bytes)\
  FOLLY_SDT(ubitrack, vision_allocate_cpu, bytes);
#endif

  /*
 * TRACEPOINT_VISION_ALLOCATE_GPU(bytes)
 * traces the allocation on GPU memory of an image
 * parameters:
 * - bytes: number of bytes that were allocated
 *
 * example:
 * TRACEPOINT_VISION_ALLOCATE_GPU(m_width*m_height*m_channels)
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_VISION_ALLOCATE_GPU(bytes)\
  if (UBITRACK_VISION_ALLOCATE_GPU_ENABLED()) {\
    UBITRACK_VISION_ALLOCATE_GPU(bytes);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_VISION_ALLOCATE_GPU(bytes)\
  ETWUbitrackAllocateGpu(bytes);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_VISION_ALLOCATE_GPU(bytes)\
  tracepoint(ubitrack, vision_allocate_gpu, bytes);
#endif


#ifdef HAVE_USDT
#define TRACEPOINT_VISION_ALLOCATE_GPU(bytes)\
  FOLLY_SDT(ubitrack, vision_allocate_gpu, bytes);
#endif

 /*
 * TRACEPOINT_VISION_GPU_UPLOAD(bytes)
 * traces the uploads to GPU memory of an image
 * parameters:
 * - bytes: number of bytes that were allocated
 *
 * example:
 * TRACEPOINT_VISION_GPU_UPLOAD(m_width*m_height*m_channels)
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_VISION_GPU_UPLOAD(bytes)\
  if (UBITRACK_VISION_GPU_UPLOAD_ENABLED()) {\
    UBITRACK_VISION_GPU_UPLOAD(bytes);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_VISION_GPU_UPLOAD(bytes)\
  ETWUbitrackGpuUpload(bytes);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_VISION_GPU_UPLOAD(bytes)\
  tracepoint(ubitrack, vision_gpu_upload, bytes);
#endif


#ifdef HAVE_USDT
#define TRACEPOINT_VISION_GPU_UPLOAD(bytes)\
  FOLLY_SDT(ubitrack, vision_gpu_upload, bytes);
#endif

 /*
 * TRACEPOINT_VISION_GPU_DOWNLOAD(bytes)
 * traces the downloas to CPU memory of an image
 * parameters:
 * - bytes: number of bytes that were allocated
 *
 * example:
 * TRACEPOINT_VISION_GPU_DOWNLOAD(m_width*m_height*m_channels)
 */
#if defined(HAVE_DTRACE) && !defined(DISABLE_DTRACE)
#define TRACEPOINT_VISION_GPU_DOWNLOAD(bytes)\
  if (UBITRACK_VISION_GPU_DOWNLOAD_ENABLED()) {\
    UBITRACK_VISION_GPU_DOWNLOAD(bytes);\
  }
#endif

#ifdef HAVE_ETW
#define TRACEPOINT_VISION_GPU_DOWNLOAD(bytes)\
  ETWUbitrackGpuDownload(bytes);
#endif

#ifdef HAVE_LTTNGUST
#define TRACEPOINT_VISION_GPU_DOWNLOAD(bytes)\
  tracepoint(ubitrack, vision_gpu_download, bytes);
#endif

#ifdef HAVE_USDT
#define TRACEPOINT_VISION_GPU_DOWNLOAD(bytes)\
  FOLLY_SDT(ubitrack, vision_gpu_download, bytes);
#endif

#else // ENABLE_EVENT_TRACING

#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_BEGIN(event_domain, event_priority, component_name, component_port)
#define TRACEPOINT_BLOCK_EVENTQUEUE_DISPATCH_END(event_domain, event_priority, component_name, component_port)
#define TRACEPOINT_MEASUREMENT_CREATE(event_domain, event_priority, component_name, component_port)
#define TRACEPOINT_MEASUREMENT_RECEIVE(event_domain, event_priority, component_name, component_port)

#define TRACEPOINT_VISION_ALLOCATE_CPU(bytes)
#define TRACEPOINT_VISION_ALLOCATE_GPU(bytes)
#define TRACEPOINT_VISION_GPU_UPLOAD(bytes)
#define TRACEPOINT_VISION_GPU_DOWNLOAD(bytes)

#endif // ENABLE_EVENT_TRACING

#endif //UBITRACK_TRACINGPROVIDER_H

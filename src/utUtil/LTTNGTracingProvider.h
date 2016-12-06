//
// Created by jack on 10.07.16.
//

#undef TRACEPOINT_PROVIDER
#define TRACEPOINT_PROVIDER ubitrack

#undef TRACEPOINT_INCLUDE
#define TRACEPOINT_INCLUDE "utUtil/LTTNGTracingProvider.h"

#if !defined(UBITRACK_LTTNGTRACINGPROVIDER_H) || defined(TRACEPOINT_HEADER_MULTI_READ)
#define UBITRACK_LTTNGTRACINGPROVIDER_H

#include <lttng/tracepoint.h>

TRACEPOINT_EVENT(
        ubitrack,
        eventqueue_dispatch_begin,
        TP_ARGS(
            unsigned int, domain,
            unsigned long long int, priority,
            const char*, component,
            const char*, port
        ),
        TP_FIELDS(
            ctf_integer(unsigned int, domain_field, domain)
            ctf_integer(unsigned long long int, priority_field, priority)
            ctf_string(component_field, component)
            ctf_string(port_field, port)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        eventqueue_dispatch_end,
        TP_ARGS(
            unsigned int, domain,
            unsigned long long int, priority,
            const char*, component,
            const char*, port
        ),
        TP_FIELDS(
            ctf_integer(unsigned int, domain_field, domain)
            ctf_integer(unsigned long long int, priority_field, priority)
            ctf_string(component_field, component)
            ctf_string(port_field, port)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        eventqueue_dispatch_discard,
        TP_ARGS(
            unsigned int, domain,
            unsigned long long int, priority,
            const char*, component,
            const char*, port
        ),
        TP_FIELDS(
            ctf_integer(unsigned int, domain_field, domain)
            ctf_integer(unsigned long long int, priority_field, priority)
            ctf_string(component_field, component)
            ctf_string(port_field, port)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        measurement_create,
        TP_ARGS(
            unsigned int, domain,
            unsigned long long int, timestamp,
            const char*, component,
            const char*, port
        ),
        TP_FIELDS(
            ctf_integer(unsigned int, domain_field, domain)
            ctf_integer(unsigned long long int, timestamp_field, timestamp)
            ctf_string(component_field, component)
            ctf_string(port_field, port)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        measurement_receive,
        TP_ARGS(
            unsigned int, domain,
            unsigned long long int, timestamp,
            const char*, component,
            const char*, port
        ),
        TP_FIELDS(
            ctf_integer(unsigned int, domain_field, domain)
            ctf_integer(unsigned long long int, timestamp_field, timestamp)
            ctf_string(component_field, component)
            ctf_string(port_field, port)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        vision_allocate_cpu,
        TP_ARGS(
        long int, size
        ),
        TP_FIELDS(
            ctf_integer(long int, size_field, size)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        vision_allocate_gpu,
        TP_ARGS(
            long int, size
        ),
        TP_FIELDS(
            ctf_integer(long int, size_field, size)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        vision_gpu_upload,
        TP_ARGS(
            long int, size
        ),
        TP_FIELDS(
            ctf_integer(long int, size_field, size)
        )
)

TRACEPOINT_EVENT(
        ubitrack,
        vision_gpu_download,
        TP_ARGS(
            long int, size
        ),
        TP_FIELDS(
            ctf_integer(long int, size_field, size)
        )
)

#endif /* UBITRACK_LTTNGTRACINGPROVIDER_H */

#include <lttng/tracepoint-event.h>

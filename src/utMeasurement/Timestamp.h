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
 * Timestamp type definition and helper functions
 * @author Florian Echtler <echtler@in.tum.de>
 */ 


#ifndef _Ubitrack_Measurement_Timestamp_INCLUDED_
#define _Ubitrack_Measurement_Timestamp_INCLUDED_

#include <string>
#include <utCore.h>

namespace Ubitrack { namespace Measurement {

/// Timestamp: nanoseconds since epoch (UNIX birth)
typedef unsigned long long int Timestamp;

/// retrieve the current system time as Timestamp
UBITRACK_EXPORT Timestamp now();

/// convert a Timestamp to a string (returns something like "Fri Mar 02 11:41:41 2007 UTC")
UBITRACK_EXPORT std::string timestampToString( Timestamp );

/// convert a Timestamp to a shorter string (returns something like "11:41:41.521021")
UBITRACK_EXPORT std::string timestampToShortString( Timestamp );

} } // namespace Ubitrack::Measurement

#endif // _Ubitrack_Measurement_Timestamp_INCLUDED_


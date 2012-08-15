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
 * Provides Ubitrack-specific functions related to the logging framework.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_UTIL_LOGGING_H_INCLUDED__
#define __UBITRACK_UTIL_LOGGING_H_INCLUDED__

#include <utCore.h>

namespace Ubitrack { namespace Util {

/** Initializes the logger */
void UBITRACK_EXPORT initLogging( const char* sConfigFile = "log4cpp.conf" );

} } // namespace Ubitrack::Util

#endif


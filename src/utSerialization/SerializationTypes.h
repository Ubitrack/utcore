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
 * @ingroup serialization
 * @file
 * Base Clases for Serialization copied/inspired by ROS Comm
 * @author Ulrich Eck <ulrich.eck@tum.de>
 */

#include "utCore.h"

// @todo check if really needed ??

#ifndef UBITRACK_SERIALIZATIONTYPES_H
#define UBITRACK_SERIALIZATIONTYPES_H

#if defined(_MSC_VER) && (_MSC_VER<1600) // MS express/studio 2008 or earlier
typedef          __int64  int64_t;
typedef unsigned __int64 uint64_t;
typedef          __int32  int32_t;
typedef unsigned __int32 uint32_t;
typedef          __int16  int16_t;
typedef unsigned __int16 uint16_t;
typedef          __int8    int8_t;
typedef unsigned __int8   uint8_t;
#else
#include <stdint.h>
#endif

#endif //UBITRACK_SERIALIZATIONTYPES_H

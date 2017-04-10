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
 * Base Clases for Serialization
 * @author Ulrich Eck <ulrich.eck@tum.de>
 */

#ifndef UBITRACK_EXCEPTION_H
#define UBITRACK_EXCEPTION_H

#include "utCore.h"
#include "utUtil/Exception.h"

#include <stdexcept>
#include <string>

#ifdef _MSC_VER
#pragma warning( disable : 4275 )
#endif

namespace Ubitrack {
namespace Serialization {

class UBITRACK_EXPORT StreamOverrunException: public Ubitrack::Util::Exception {
public:
    StreamOverrunException(const std::string& what) : Exception(what) {}
};

UBITRACK_EXPORT void throwStreamOverrun();

} // Serialization
} // Ubitrack

#endif //UBITRACK_EXCEPTION_H

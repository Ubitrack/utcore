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


#ifndef UBITRACK_SERIALIZATIONFORMAT_H
#define UBITRACK_SERIALIZATIONFORMAT_H

#include "utSerialization/Exception.h"
#include "utSerialization/SerializationTypes.h"
#include "utSerialization/SerializationTraits.h"

#include <boost/array.hpp>
#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>

namespace Ubitrack {
namespace Serialization {
namespace Format {

/**
 * \brief Templated serialization class.
 *
 * Specializing the SerializationFormat class is the only thing you need to do to get the serialization system
 * to work with a type.
 */
template<typename T>
struct SerializationFormat {
  /**
   * \brief Write an object to the stream.
   */
  template<typename Stream>
  inline static void write(Stream& stream, typename boost::call_traits<T>::param_type t)
  {
  }

  /**
   * \brief Read an object from the stream.
   */
  template<typename Stream>
  inline static void read(Stream& stream, typename boost::call_traits<T>::reference t)
  {
  }

  /**
   * \brief Determine the maximum serialized length of an object.
   */
  inline static uint32_t maxSerializedLength(typename boost::call_traits<T>::param_type t)
  {
      return 0;
  }
};

} // Format
} // Serialization
} // Ubitrack

#endif //UBITRACK_SERIALIZATIONFORMAT_H

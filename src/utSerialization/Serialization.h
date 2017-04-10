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


#ifndef UBITRACK_SERIALIZATION_H
#define UBITRACK_SERIALIZATION_H

#include "utSerialization/Exception.h"
#include "utSerialization/SerializationTypes.h"
#include "utSerialization/SerializationTraits.h"
#include <stdexcept>
#include <string>


#include <boost/array.hpp>
#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>



namespace Ubitrack {
namespace Serialization {

/**
 * \brief Templated serialization class.
 *
 * Base Serializer Class - needs a format to serialize/deserialize Objects of type T to a stream of type Stream
 */
template<typename T, typename F>
struct Serializer {

  /**
   * \brief Write an object to the stream.
   */
  template<typename Stream>
  inline static void write(Stream& stream, typename boost::call_traits<T>::param_type t)
  {
      F::write(stream, t);
  }

  /**
   * \brief Read an object from the stream.
   */
  template<typename Stream>
  inline static void read(Stream& stream, typename boost::call_traits<T>::reference t)
  {
      F::read(stream, t);
  }

  /**
   * \brief Determine the maximum serialized length of an object.
   */
  inline static uint32_t maxSerializedLength(typename boost::call_traits<T>::param_type t)
  {
      return F::maxSerializedLength(t);
  }
};

/**
 * \brief Serialize an object.
 */
template<typename T, typename F, typename Stream>
inline void serialize(Stream& stream, const T& t)
{
    Serializer<T, F>::write(stream, t);
}

/**
 * \brief Deserialize an object.
 */
template<typename T, typename F, typename Stream>
inline void deserialize(Stream& stream, T& t)
{
    Serializer<T, F>::read(stream, t);
}

/**
 * \brief Determine the maximum serialized length of an object
 */
template<typename T, typename F>
inline uint32_t maxSerializationLength(const T& t)
{
    return Serializer<T, F>::maxSerializedLength(t);
}


} // Serialization
} // Ubitrack

#endif //UBITRACK_SERIALIZATION_H

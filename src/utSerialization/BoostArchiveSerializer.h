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
 * Boost-Binary Serialization
 * @author Ulrich Eck <ulrich.eck@tum.de>
 */

#ifndef UBITRACK_BOOSTBINARYSERIALIZER_H
#define UBITRACK_BOOSTBINARYSERIALIZER_H


#include "utSerialization/Serialization.h"
#include "utSerialization/SerializationFormat.h"

#include "utMeasurement/Measurement.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include <boost/array.hpp>
#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace Ubitrack {
namespace Serialization {
namespace BoostArchive {


template<typename T>
struct BoostArchiveSerializationFormat
{
  template<typename Stream>
  inline static void write(Stream& stream,  typename boost::call_traits<T>::param_type t)
  {
      stream << t;
  }

  template<typename Stream>
  inline static void read(Stream& stream, typename boost::call_traits<T>::reference t)
  {
      stream >> t;
  }

  inline static uint32_t maxSerializedLength(typename boost::call_traits<T>::param_type t)
  {
      // only way to find out is to either specialize for all types
      // or to serialize into a stringstream  to know for shure ...
      return 0; // Not Supported for now
  }
};


/**
 * \brief Serialize an object.  Stream here should normally be a boost::archive::binary_oarchive
 */
template<typename T, typename Stream>
inline void serialize(Stream& stream, const T& t)
{
    Serializer<T, BoostArchiveSerializationFormat<T> >::write(stream, t);
}

/**
 * \brief Deserialize an object.  Stream here should normally be a boost::archive::binary_iarchive
 */
template<typename T, typename Stream>
inline void deserialize(Stream& stream, T& t)
{
    Serializer<T, BoostArchiveSerializationFormat<T> >::read(stream, t);
}

/**
 * \brief Determine the serialized length of an object
 */
template<typename T>
inline uint32_t maxSerializationLength(const T& t)
{
    return Serializer<T, BoostArchiveSerializationFormat<T> >::maxSerializedLength(t);
}



} // BoostArchive
} // Serialization
} // Ubitrack

#endif //UBITRACK_BOOSTBINARYSERIALIZER_H

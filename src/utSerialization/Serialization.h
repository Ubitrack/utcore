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
 * Frontend for Serialization of Ubitrack Types
 * @author Ulrich Eck <ulrich.eck@tum.de>
 */


#ifndef UBITRACK_SERIALIZATION_H
#define UBITRACK_SERIALIZATION_H

#include "utSerialization/BoostArchiveSerializer.h"
#include "utSerialization/MsgpackSerializer.h"


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


enum SerializationProtocol {
  PROTOCOL_UNKNOWN = 0,
  PROTOCOL_BOOST_TEXT,
  PROTOCOL_BOOST_BINARY,
  PROTOCOL_MSGPACK
};


/**
 * \brief Serialize an object.
 */
template<typename T, typename Stream>
inline void serialize(const SerializationProtocol p, Stream& stream, const T& t)
{
    switch(p) {
    case PROTOCOL_BOOST_TEXT:
    {
        boost::archive::text_oarchive out_archive_t( stream );
        BoostArchive::serialize(out_archive_t, t);
    }
        break;
    case PROTOCOL_BOOST_BINARY:
    {
        boost::archive::binary_oarchive out_archive_b( stream );
        BoostArchive::serialize(out_archive_b, t);
    }
        break;
    case PROTOCOL_MSGPACK:
        MsgpackArchive::serialize(stream, t);
        break;
    default:
        UBITRACK_THROW("Unknown Serialization Protocol");
    }
}

/**
 * \brief Deserialize an object using a protocol (format).
 */
template<typename T, typename Stream>
inline void deserialize(const SerializationProtocol p, Stream& stream, T& t)
{
    switch(p) {
    case PROTOCOL_BOOST_TEXT:
    {
        boost::archive::text_iarchive in_archive_t( stream );
        BoostArchive::deserialize(in_archive_t, t);
    }
        break;
    case PROTOCOL_BOOST_BINARY:
    {
        boost::archive::binary_iarchive in_archive_b( stream );
        BoostArchive::deserialize(in_archive_b, t);
    }
        break;
    case PROTOCOL_MSGPACK:
        MsgpackArchive::deserialize(stream, t);
        break;
    default:
        UBITRACK_THROW("Unknown Serialization Protocol");
    }
}


} // Serialization
} // Ubitrack

#endif //UBITRACK_SERIALIZATION_H

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

#ifndef UBITRACK_SERIALIZATIONTRAITS_H
#define UBITRACK_SERIALIZATIONTRAITS_H


#include <string>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include "utSerialization/SerializationTypes.h"

namespace Ubitrack {
namespace Serialization {
namespace Traits {
/**
 * \brief Base type for compile-type true/false tests.  Compatible with Boost.MPL.  classes inheriting from this type
 * are \b true values.
 */
struct TrueType {
  static const bool value = true;
  typedef TrueType type;
};

/**
 * \brief Base type for compile-type true/false tests.  Compatible with Boost.MPL.  classes inheriting from this type
 * are \b false values.
 */
struct FalseType {
  static const bool value = false;
  typedef FalseType type;
};

/**
 * \brief A simple datatype is one that can be memcpy'd directly in array form, i.e. it's a POD, fixed-size type and
 * sizeof(M) == sum(serializationLength(M:a...))
 */
template<typename M>
struct IsSimple: public FalseType {};
/**
 * \brief A fixed-size datatype is one whose size is constant, i.e. it has no variable-length arrays or strings
 */
template<typename M>
struct IsFixedSize: public FalseType {};
/**
 * \brief HasHeader informs whether or not there is a header that gets serialized as the first thing in the message
 */
template<typename M>
struct HasHeader: public FalseType {};

/**
 * \brief Am I message or not
 */
template<typename M>
struct IsMessage: public FalseType {};

/**
 * \brief returns IsSimple<M>::value;
 */
template<typename M>
inline bool isSimple()
{
    return IsSimple<typename boost::remove_reference<typename boost::remove_const<M>::type>::type>::value;
}

/**
 * \brief returns IsFixedSize<M>::value;
 */
template<typename M>
inline bool isFixedSize()
{
    return IsFixedSize<typename boost::remove_reference<typename boost::remove_const<M>::type>::type>::value;
}

#define UBITRACK_CREATE_SIMPLE_TRAITS(Type) \
    template<> struct IsSimple<Type> : public TrueType {}; \
    template<> struct IsFixedSize<Type> : public TrueType {};

UBITRACK_CREATE_SIMPLE_TRAITS(uint8_t);
UBITRACK_CREATE_SIMPLE_TRAITS(int8_t);
UBITRACK_CREATE_SIMPLE_TRAITS(uint16_t);
UBITRACK_CREATE_SIMPLE_TRAITS(int16_t);
UBITRACK_CREATE_SIMPLE_TRAITS(uint32_t);
UBITRACK_CREATE_SIMPLE_TRAITS(int32_t);
UBITRACK_CREATE_SIMPLE_TRAITS(uint64_t);
UBITRACK_CREATE_SIMPLE_TRAITS(int64_t);
UBITRACK_CREATE_SIMPLE_TRAITS(float);
UBITRACK_CREATE_SIMPLE_TRAITS(double);

// because std::vector<bool> is not a true vector, bool is not a simple type
template<>
struct IsFixedSize<bool>: public TrueType {};

} // Traits
} // Serialization
} // Ubitrack

#endif //UBITRACK_SERIALIZATIONTRAITS_H

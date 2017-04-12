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
 * Custom (probably ROS Compatible) Serialization
 * @author Ulrich Eck <ulrich.eck@tum.de>
 */


#ifndef UBITRACK_BINARYSERIALIZATION_H
#define UBITRACK_BINARYSERIALIZATION_H

#include "utSerialization/Serialization.h.h"
#include "utSerialization/SerializationFormat.h"

#include <boost/array.hpp>
#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>

#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>

namespace Ubitrack {
namespace Serialization {
namespace ROSBinary {

namespace mpl = boost::mpl;
namespace mt = Ubitrack::Serialization::Traits;

/**
 * \brief Templated serialization class.  Default implementation provides backwards compatibility with
 * old message types.
 *
 * Specializing the Serializer class is the only thing you need to do to get the ROS serialization system
 * to work with a type.
 */
template<typename T>
struct ROSBinarySerializationFormat {
  /**
   * \brief Write an object to the stream.  Normally the stream passed in here will be a ros::serialization::OStream
   */
  template<typename Stream>
  inline static void write(Stream& stream, typename boost::call_traits<T>::param_type t)
  {
      t.serialize(stream.getData(), 0);
  }

  /**
   * \brief Read an object from the stream.  Normally the stream passed in here will be a ros::serialization::IStream
   */
  template<typename Stream>
  inline static void read(Stream& stream, typename boost::call_traits<T>::reference t)
  {
      t.deserialize(stream.getData());
  }

  /**
   * \brief Determine the serialized length of an object.
   */
  inline static uint32_t maxSerializedLength(typename boost::call_traits<T>::param_type t)
  {
      return t.maxSerializationLength();
  }
};

/**
 * \brief Serialize an object.  Stream here should normally be a ros::serialization::OStream
 */
template<typename T, typename Stream>
inline void serialize(Stream& stream, const T& t)
{
    Serializer<T, ROSBinarySerializationFormat<T> >::write(stream, t);
}

/**
 * \brief Deserialize an object.  Stream here should normally be a ros::serialization::IStream
 */
template<typename T, typename Stream>
inline void deserialize(Stream& stream, T& t)
{
    Serializer<T, ROSBinarySerializationFormat<T> >::read(stream, t);
}

/**
 * \brief Determine the serialized length of an object
 */
template<typename T>
inline uint32_t maxSerializationLength(const T& t)
{
    return Serializer<T, ROSBinarySerializationFormat<T> >::maxSerializedLength(t);
}

#define ROSBINARY_CREATE_SIMPLE_SERIALIZER(Type) \
  template<> struct ROSBinarySerializationFormat<Type> \
  { \
    template<typename Stream> inline static void write(Stream& stream, const Type v) \
    { \
      *reinterpret_cast<Type*>(stream.advance(sizeof(v))) = v; \
    } \
    \
    template<typename Stream> inline static void read(Stream& stream, Type& v) \
    { \
      v = *reinterpret_cast<Type*>(stream.advance(sizeof(v))); \
    } \
    \
    inline static uint32_t maxSerializedLength(const Type&) \
    { \
      return sizeof(Type); \
    } \
  };

ROSBINARY_CREATE_SIMPLE_SERIALIZER(uint8_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(int8_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(uint16_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(int16_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(uint32_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(int32_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(uint64_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(int64_t);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(float);
ROSBINARY_CREATE_SIMPLE_SERIALIZER(double);

/**
 * \brief Serializer specialized for bool (serialized as uint8)
 */
template<>
struct ROSBinarySerializationFormat<bool> {
  template<typename Stream>
  inline static void write(Stream& stream, const bool v)
  {
      uint8_t b = (uint8_t) v;
      *reinterpret_cast<uint8_t*>(stream.advance(1)) = b;
  }

  template<typename Stream>
  inline static void read(Stream& stream, bool& v)
  {
      uint8_t b;
      b = *reinterpret_cast<uint8_t*>(stream.advance(1));
      v = (bool) b;
  }

  inline static uint32_t maxSerializedLength(bool)
  {
      return 1;
  }
};

/**
 * \brief  Serializer specialized for std::string
 */
template<class ContainerAllocator>
struct ROSBinarySerializationFormat<std::basic_string<char, std::char_traits<char>, ContainerAllocator> > {
  typedef std::basic_string<char, std::char_traits<char>, ContainerAllocator> StringType;

  template<typename Stream>
  inline static void write(Stream& stream, const StringType& str)
  {
      size_t len = str.size();
      stream.next((uint32_t) len);

      if (len>0) {
          memcpy(stream.advance((uint32_t) len), str.data(), len);
      }
  }

  template<typename Stream>
  inline static void read(Stream& stream, StringType& str)
  {
      uint32_t len;
      stream.next(len);
      if (len>0) {
          str = StringType((char*) stream.advance(len), len);
      }
      else {
          str.clear();
      }
  }

  inline static uint32_t maxSerializedLength(const StringType& str)
  {
      return 4+(uint32_t) str.size();
  }
};

/**
 * \brief Vector serializer.  Default implementation does nothing
 */
template<typename T, class ContainerAllocator, class Enabled = void>
struct VectorSerializer {};

/**
 * \brief Vector serializer, specialized for non-fixed-size, non-simple types
 */
template<typename T, class ContainerAllocator>
struct VectorSerializer<T, ContainerAllocator, typename boost::disable_if<mt::IsFixedSize<T> >::type> {
  typedef std::vector<T, typename ContainerAllocator::template rebind<T>::other> VecType;
  typedef typename VecType::iterator IteratorType;
  typedef typename VecType::const_iterator ConstIteratorType;

  template<typename Stream>
  inline static void write(Stream& stream, const VecType& v)
  {
      stream.next((uint32_t) v.size());
      ConstIteratorType it = v.begin();
      ConstIteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  template<typename Stream>
  inline static void read(Stream& stream, VecType& v)
  {
      uint32_t len;
      stream.next(len);
      v.resize(len);
      IteratorType it = v.begin();
      IteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  inline static uint32_t maxSerializedLength(const VecType& v)
  {
      uint32_t size = 4;
      ConstIteratorType it = v.begin();
      ConstIteratorType end = v.end();
      for (; it!=end; ++it) {
          size += maxSerializationLength(*it);
      }

      return size;
  }
};

/**
 * \brief Vector serializer, specialized for fixed-size simple types
 */
template<typename T, class ContainerAllocator>
struct VectorSerializer<T, ContainerAllocator, typename boost::enable_if<mt::IsSimple<T> >::type> {
  typedef std::vector<T, typename ContainerAllocator::template rebind<T>::other> VecType;
  typedef typename VecType::iterator IteratorType;
  typedef typename VecType::const_iterator ConstIteratorType;

  template<typename Stream>
  inline static void write(Stream& stream, const VecType& v)
  {
      uint32_t len = (uint32_t) v.size();
      stream.next(len);
      if (!v.empty()) {
          const uint32_t data_len = len*(uint32_t) sizeof(T);
          memcpy(stream.advance(data_len), &v.front(), data_len);
      }
  }

  template<typename Stream>
  inline static void read(Stream& stream, VecType& v)
  {
      uint32_t len;
      stream.next(len);
      v.resize(len);

      if (len>0) {
          const uint32_t data_len = (uint32_t) sizeof(T)*len;
          memcpy(&v.front(), stream.advance(data_len), data_len);
      }
  }

  inline static uint32_t maxSerializedLength(const VecType& v)
  {
      return 4+v.size()*(uint32_t) sizeof(T);
  }
};

/**
 * \brief Vector serializer, specialized for fixed-size non-simple types
 */
template<typename T, class ContainerAllocator>
struct VectorSerializer<T,
                        ContainerAllocator,
                        typename boost::enable_if<mpl::and_<mt::IsFixedSize<T>,
                                                            mpl::not_<mt::IsSimple<T> > > >::type> {
  typedef std::vector<T, typename ContainerAllocator::template rebind<T>::other> VecType;
  typedef typename VecType::iterator IteratorType;
  typedef typename VecType::const_iterator ConstIteratorType;

  template<typename Stream>
  inline static void write(Stream& stream, const VecType& v)
  {
      stream.next((uint32_t) v.size());
      ConstIteratorType it = v.begin();
      ConstIteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  template<typename Stream>
  inline static void read(Stream& stream, VecType& v)
  {
      uint32_t len;
      stream.next(len);
      v.resize(len);
      IteratorType it = v.begin();
      IteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  inline static uint32_t maxSerializedLength(const VecType& v)
  {
      uint32_t size = 4;
      if (!v.empty()) {
          uint32_t len_each = maxSerializationLength(v.front());
          size += len_each*(uint32_t) v.size();
      }

      return size;
  }
};

/**
 * \brief serialize version for std::vector
 */
template<typename T, class ContainerAllocator, typename Stream>
inline void serialize(Stream& stream, const std::vector<T, ContainerAllocator>& t)
{
    VectorSerializer<T, ContainerAllocator>::write(stream, t);
}

/**
 * \brief deserialize version for std::vector
 */
template<typename T, class ContainerAllocator, typename Stream>
inline void deserialize(Stream& stream, std::vector<T, ContainerAllocator>& t)
{
    VectorSerializer<T, ContainerAllocator>::read(stream, t);
}

/**
 * \brief maxSerializationLength version for std::vector
 */
template<typename T, class ContainerAllocator>
inline uint32_t maxSerializationLength(const std::vector<T, ContainerAllocator>& t)
{
    return VectorSerializer<T, ContainerAllocator>::serializedLength(t);
}

/**
 * \brief Array serializer, default implementation does nothing
 */
template<typename T, size_t N, class Enabled = void>
struct ArraySerializer {};

/**
 * \brief Array serializer, specialized for non-fixed-size, non-simple types
 */
template<typename T, size_t N>
struct ArraySerializer<T, N, typename boost::disable_if<mt::IsFixedSize<T> >::type> {
  typedef boost::array<T, N> ArrayType;
  typedef typename ArrayType::iterator IteratorType;
  typedef typename ArrayType::const_iterator ConstIteratorType;

  template<typename Stream>
  inline static void write(Stream& stream, const ArrayType& v)
  {
      ConstIteratorType it = v.begin();
      ConstIteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  template<typename Stream>
  inline static void read(Stream& stream, ArrayType& v)
  {
      IteratorType it = v.begin();
      IteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  inline static uint32_t maxSerializedLength(const ArrayType& v)
  {
      uint32_t size = 0;
      ConstIteratorType it = v.begin();
      ConstIteratorType end = v.end();
      for (; it!=end; ++it) {
          size += maxSerializationLength(*it);
      }

      return size;
  }
};

/**
 * \brief Array serializer, specialized for fixed-size, simple types
 */
template<typename T, size_t N>
struct ArraySerializer<T, N, typename boost::enable_if<mt::IsSimple<T> >::type> {
  typedef boost::array<T, N> ArrayType;
  typedef typename ArrayType::iterator IteratorType;
  typedef typename ArrayType::const_iterator ConstIteratorType;

  template<typename Stream>
  inline static void write(Stream& stream, const ArrayType& v)
  {
      const uint32_t data_len = N*sizeof(T);
      memcpy(stream.advance(data_len), &v.front(), data_len);
  }

  template<typename Stream>
  inline static void read(Stream& stream, ArrayType& v)
  {
      const uint32_t data_len = N*sizeof(T);
      memcpy(&v.front(), stream.advance(data_len), data_len);
  }

  inline static uint32_t maxSerializedLength(const ArrayType&)
  {
      return N*sizeof(T);
  }
};

/**
 * \brief Array serializer, specialized for fixed-size, non-simple types
 */
template<typename T, size_t N>
struct ArraySerializer<T,
                       N,
                       typename boost::enable_if<mpl::and_<mt::IsFixedSize<T>,
                                                           mpl::not_<mt::IsSimple<T> > > >::type> {
  typedef boost::array<T, N> ArrayType;
  typedef typename ArrayType::iterator IteratorType;
  typedef typename ArrayType::const_iterator ConstIteratorType;

  template<typename Stream>
  inline static void write(Stream& stream, const ArrayType& v)
  {
      ConstIteratorType it = v.begin();
      ConstIteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  template<typename Stream>
  inline static void read(Stream& stream, ArrayType& v)
  {
      IteratorType it = v.begin();
      IteratorType end = v.end();
      for (; it!=end; ++it) {
          stream.next(*it);
      }
  }

  inline static uint32_t maxSerializedLength(const ArrayType& v)
  {
      return maxSerializationLength(v.front())*N;
  }
};

/**
 * \brief serialize version for boost::array
 */
template<typename T, size_t N, typename Stream>
inline void serialize(Stream& stream, const boost::array<T, N>& t)
{
    ArraySerializer<T, N>::write(stream, t);
}

/**
 * \brief deserialize version for boost::array
 */
template<typename T, size_t N, typename Stream>
inline void deserialize(Stream& stream, boost::array<T, N>& t)
{
    ArraySerializer<T, N>::read(stream, t);
}

/**
 * \brief maxSerializationLength version for boost::array
 */
template<typename T, size_t N>
inline uint32_t maxSerializationLength(const boost::array<T, N>& t)
{
    return ArraySerializer<T, N>::serializedLength(t);
}


/**
 * \brief Enum
 */
namespace stream_types {
enum StreamType {
  Input,
  Output,
  Length
};
}
typedef stream_types::StreamType StreamType;

/**
 * \brief Stream base-class, provides common functionality for IStream and OStream
 */
struct UBITRACK_EXPORT Stream {
    /*
     * \brief Returns a pointer to the current position of the stream
     */
    inline uint8_t* getData() { return data_; }
    /**
     * \brief Advances the stream, checking bounds, and returns a pointer to the position before it
     * was advanced.
     * \throws StreamOverrunException if len would take this stream past the end of its buffer
     */
    ROSBINARY_FORCE_INLINE uint8_t
    *
    advance(uint32_t
    len)
    {
        uint8_t*old_data = data_;
        data_ += len;
        if (data_>end_) {
            // Throwing directly here causes a significant speed hit due to the extra code generated
            // for the throw statement
            throwStreamOverrun();
        }
        return old_data;
    }

    /**
     * \brief Returns the amount of space left in the stream
     */
    inline uint32_t getLength() { return (uint32_t) (end_-data_); }

protected:
    Stream(uint8_t* _data, uint32_t _count)
            :data_(_data), end_(_data+_count) {}

private:
    uint8_t* data_;
    uint8_t* end_;
};

/**
 * \brief Input stream
 */
struct UBITRACK_EXPORT IStream
        : public Stream {
  static const StreamType stream_type = stream_types::Input;

  IStream(uint8_t* data, uint32_t count)
          :Stream(data, count)
  {
  }

/**
 * \brief Deserialize an item from this input stream
 */
  template<typename T>
  inline void next(T& t)
  {
      deserialize(*this, t);
  }

  template<typename T>
  inline IStream& operator>>(T& t)
  {
      deserialize(*this, t);
      return *this;
  }
};

/**
 * \brief Output stream
 */
struct UBITRACK_EXPORT OStream
        : public Stream {
  static const StreamType stream_type = stream_types::Output;

  OStream(uint8_t* data, uint32_t count)
          :Stream(data, count)
  {
  }

/**
 * \brief Serialize an item to this output stream
 */
  template<typename T>
  inline void next(const T& t)
  {
      serialize(*this, t);
  }

  template<typename T>
  inline OStream& operator<<(const T& t)
  {
      serialize(*this, t);
      return *this;
  }
};

/**
 * \brief Length stream
 *
 * LStream is not what you would normally think of as a stream, but it is used in order to support
 * allinone serializers.
 */
struct UBITRACK_EXPORT LStream {
    static const StreamType stream_type = stream_types::Length;

    LStream()
            :count_(0) {}

    /**
     * \brief Add the length of an item to this length stream
     */
    template<typename T>
    inline void next(const T& t)
    {
        count_ += maxSerializationLength(t);
    }

    /**
     * \brief increment the length by len
     */
    inline uint32_t advance(uint32_t len)
    {
        uint32_t old = count_;
        count_ += len;
        return old;
    }

    /**
     * \brief Get the total length of this stream
     */
    inline uint32_t getLength() { return count_; }

private:
    uint32_t count_;
};

} // namespace ROSBinary
} // namespace Serialization
} // Ubitrack

#endif //UBITRACK_BINARYSERIALIZATION_H


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
 * Msgpack Serialization
 * @author Ulrich Eck <ulrich.eck@tum.de>
 */

// good tutorial: https://github.com/msgpack/msgpack-c/wiki/v2_0_cpp_adaptor


#ifndef UBITRACK_MSGPACKSERIALIZER_H
#define UBITRACK_MSGPACKSERIALIZER_H

#include "utSerialization/BaseSerializer.h"
#include "utSerialization/SerializationFormat.h"

#include "utMeasurement/Measurement.h"

#include <boost/array.hpp>
#include <boost/call_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <msgpack.hpp>

namespace Ubitrack {
namespace Serialization {
namespace MsgpackArchive {

template<typename T>
struct MsgpackSerializationFormat {
  template<typename Stream>
  inline static void write(Stream& stream, typename boost::call_traits<T>::param_type t)
  {
      msgpack::pack(stream, t);
  }

  template<typename Stream>
  inline static void write(msgpack::packer<Stream>& pac, typename boost::call_traits<T>::param_type t)
  {
      pac.pack(t);
  }

  template<typename Stream>
  inline static void read(Stream& pac, typename boost::call_traits<T>::reference t)
  {
      msgpack::object_handle oh;
      if (pac.next(oh)) {
          msgpack::object obj = oh.get();
          obj.convert<T>(t);
      } else {
          // throw ??
      }
  }


  inline static uint32_t maxSerializedLength(typename boost::call_traits<T>::param_type t)
  {
      return 0;
  }
};


/**
 * \brief Serialize an object.  Stream here should normally be a boost::archive::binary_oarchive
 */
template<typename T, typename Stream>
inline void serialize(Stream& stream, const T& t)
{
    BaseSerializer<T, MsgpackSerializationFormat<T> >::write(stream, t);
}

/**
 * \brief Deserialize an object.  Stream here should normally be a boost::archive::binary_iarchive
 */
template<typename T, typename Stream>
inline void deserialize(Stream& stream, T& t)
{
    BaseSerializer<T, MsgpackSerializationFormat<T> >::read(stream, t);
}

/**
 * \brief Determine the serialized length of an object
 */
template<typename T>
inline uint32_t maxSerializationLength(const T& t)
{
    return BaseSerializer<T, MsgpackSerializationFormat<T> >::maxSerializedLength(t);
}

} // MsgpackArchive
} // Serialization
} // Ubitrack


namespace msgpack {
MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
namespace adaptor {

/*
 * Measurements of type T
 */

#if !defined(MSGPACK_USE_CPP03)
template<typename T>
struct as<Ubitrack::Measurement::Measurement<T>, typename std::enable_if<msgpack::has_as<T>::value>::type> {
  Ubitrack::Measurement::Measurement<T> operator()(msgpack::object const& o) const
  {
      if (o.is_nil()) return nullptr;
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size!=2) throw msgpack::type_error();
      return Ubitrack::Measurement::Measurement<T>(
              o.via.array.ptr[0].as<unsigned long long>(),
              boost::shared_ptr<T>(o.via.array.ptr[0].as<T>())
      );
  }
};
#endif // !defined(MSGPACK_USE_CPP03)

template<typename T>
struct convert<Ubitrack::Measurement::Measurement<T> > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Measurement::Measurement<T>& v) const
  {
      if (o.is_nil()) {
          v.invalidate();
          v.reset();
      }
      else {
          if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size!=2) throw msgpack::type_error();
          T value;
          msgpack::adaptor::convert<T>()(o.via.array.ptr[1], value);
		  unsigned long long ts;
		  msgpack::adaptor::convert<unsigned long long>()(o.via.array.ptr[0], ts);
          v = Ubitrack::Measurement::Measurement<T>(ts, value);
      }
      return o;
  }
};

template<typename T>
struct pack<Ubitrack::Measurement::Measurement<T> > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Measurement::Measurement<T>& v) const
  {
      if (v) {
          o.pack_array(2);
          o.pack((unsigned long long) v.time());
          o.pack(*v);
      }
      else {
          o.pack_nil();
      }
      return o;
  }
};

template<typename T>
struct object_with_zone<Ubitrack::Measurement::Measurement<T> > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Measurement::Measurement<T>& v) const
  {
      if (v) {
          o.type = type::ARRAY;
          o.via.array.size = 2;
          o.via.array.ptr = static_cast<msgpack::object*>(
                  o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
          o.via.array.ptr[0] = msgpack::object((unsigned long long) v.time(), o.zone);
          o.via.array.ptr[1] = msgpack::object(*v, o.zone);
      }
      else {
          o.type = msgpack::type::NIL;
      }
  }
};

/*
 * Ubitrack::Math::Scalar<T>
 */
#if !defined(MSGPACK_USE_CPP03)
template<typename T>
struct as<Ubitrack::Math::Scalar<T>, typename std::enable_if<msgpack::has_as<T>::value>::type> {
  Ubitrack::Math::Scalar<T> operator()(msgpack::object const& o) const
  {
      return Ubitrack::Math::Scalar<T>(o.as<T>());
  }
};
#endif

template<typename T>
struct convert<Ubitrack::Math::Scalar<T> > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::Scalar<T>& v) const
  {
	  msgpack::adaptor::convert<T>()(o, v.m_value);
      return o;
  }
};

template<typename T>
struct pack<Ubitrack::Math::Scalar<T> > {
  template<typename Stream>
  packer<Stream>& operator()(msgpack::packer<Stream>& o, Ubitrack::Math::Scalar<T> const& v) const
  {
      o.pack(v.m_value);
      return o;
  }
};

template<typename T>
struct object_with_zone<Ubitrack::Math::Scalar<T> > {
  void operator()(msgpack::object::with_zone& o, Ubitrack::Math::Scalar<T> const& v) const
  {
      msgpack::adaptor::object_with_zone<T>()(o, v.m_value);
  }
};


/*
 * Ubitrack::Math::Vector<T, N>
 */
#if !defined(MSGPACK_USE_CPP03)
template<typename T, std::size_t N>
struct as<Ubitrack::Math::Vector<T, N>, typename std::enable_if<msgpack::has_as<T>::value>::type> {
  Ubitrack::Math::Vector<T, N> operator()(msgpack::object const& o) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      std::size_t num_elements = o.via.array.size;
      if (N>0) {
          if (num_elements!=N) throw msgpack::type_error();
      }
      Ubitrack::Math::Vector<T, N> result(num_elements);
      for (std::size_t i = 0; i<num_elements; ++i) {
        result(i) = o.via.array.ptr[i].as<T>();
      }
      return result;
  }
};
#endif // !defined(MSGPACK_USE_CPP03)

template<typename T, std::size_t N>
struct convert<Ubitrack::Math::Vector<T, N> > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::Vector<T, N>& v) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      std::size_t num_elements = o.via.array.size;
      if (N>0) {
          if (num_elements!=N) throw msgpack::type_error();
      }
      for (std::size_t i = 0; i<num_elements; ++i) {
		  msgpack::adaptor::convert<T>()(o.via.array.ptr[i], v(i));
          //v(i) = o.via.array.ptr[i].as<T>();
      }
      return o;
  }
};

template<typename T, std::size_t N>
struct pack<Ubitrack::Math::Vector<T, N> > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::Vector<T, N>& v) const
  {
      std::size_t num_elements = N;
      if (num_elements == 0) {
          num_elements = v.size();
      }
      o.pack_array((uint32_t)num_elements);
      for (std::size_t i = 0; i<num_elements; ++i) {
          o.pack(v(i));
      }
      return o;
  }
};

template<typename T, std::size_t N>
struct object_with_zone<Ubitrack::Math::Vector<T, N> > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::Vector<T, N>& v) const
  {
      std::size_t num_elements = N;
      if (num_elements == 0) {
          num_elements = v.size();
      }
      o.type = type::ARRAY;
      o.via.array.size = (uint32_t)num_elements;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      for (std::size_t i = 0; i<num_elements; ++i) {
          o.via.array.ptr[0] = msgpack::object(v(i), o.zone);
      }
  }
};


/*
 * Ubitrack::Math::Matrix<T, M, N>
 */
#if !defined(MSGPACK_USE_CPP03)
template<typename T, std::size_t M, std::size_t N>
struct as<Ubitrack::Math::Matrix<T, M, N>, typename std::enable_if<msgpack::has_as<T>::value>::type> {
  Ubitrack::Math::Matrix<T, M, N> operator()(msgpack::object const& o) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      std::size_t num_elements = o.via.array.size;
      if ((M>0) && (N>0)) {
          if (num_elements!=M*N) throw msgpack::type_error();
      } else {
          // cannot unpack dynamically sized matrix without knowing dimensions ...
          throw msgpack::type_error();
      }
      Ubitrack::Math::Matrix<T, M, N> result(M, N);
      for (std::size_t i = 0; i<M; ++i) {
          for (std::size_t j = 0; j<N; ++j) {
              std::size_t idx = i*M + j;
              result(i,j) = o.via.array.ptr[idx].as<T>();
          }
      }
      return result;
  }
};
#endif // !defined(MSGPACK_USE_CPP03)

template<typename T, std::size_t M, std::size_t N>
struct convert<Ubitrack::Math::Matrix<T, M, N> > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::Matrix<T, M, N>& v) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      std::size_t num_elements = o.via.array.size;
      if ((M>0) && (N>0)) {
          if (num_elements!=M*N) throw msgpack::type_error();
      } else {
          // cannot unpack dynamically sized matrix without knowing dimensions ...
          throw msgpack::type_error();
      }
      for (std::size_t i = 0; i<M; ++i) {
          for (std::size_t j = 0; j<N; ++j) {
              std::size_t idx = i*M + j;
			  msgpack::adaptor::convert<T>()(o.via.array.ptr[idx], v(i, j));
              //v(i,j) = o.via.array.ptr[idx].as<T>();
          }
      }
      return o;
  }
};

template<typename T, std::size_t M, std::size_t N>
struct pack<Ubitrack::Math::Matrix<T, M, N> > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::Matrix<T, M, N>& v) const
  {
      std::size_t num_elements = M*N;
      if (num_elements == 0) {
          // cannot pack dynamically sized matrix without also serializing dimensions ...
          throw msgpack::type_error();
      }
      o.pack_array((uint32_t)num_elements);
      for (std::size_t i = 0; i<M; ++i) {
          for (std::size_t j = 0; j<N; ++j) {
              std::size_t idx = i*M + j;
              o.pack(v(i,j));
          }
      }
      return o;
  }
};

template<typename T, std::size_t M, std::size_t N>
struct object_with_zone<Ubitrack::Math::Matrix<T, M, N> > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::Matrix<T, M, N>& v) const
  {
      std::size_t num_elements = M*N;
      if (num_elements == 0) {
          // cannot pack dynamically sized matrix without also serializing dimensions ...
          throw msgpack::type_error();
      }
      o.type = type::ARRAY;
      o.via.array.size = (uint32_t)num_elements;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      for (std::size_t i = 0; i<M; ++i) {
          for (std::size_t j = 0; j<N; ++j) {
              std::size_t idx = i*M + j;
              o.via.array.ptr[idx] = msgpack::object(v(i,j), o.zone);
          }
      }
  }
};





/*
 * Ubitrack::Math::Quaternion
 */
#if !defined(MSGPACK_USE_CPP03)
template<>
struct as<Ubitrack::Math::Quaternion> {
  Ubitrack::Math::Quaternion operator()(msgpack::object const& o) const
  {
      return Ubitrack::Math::Quaternion::fromVector(o.as<Ubitrack::Math::Vector<double, 4> >());
  }
};
#endif 

template<>
struct convert<Ubitrack::Math::Quaternion > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::Quaternion& v) const
  {
	  Ubitrack::Math::Vector<double, 4> vec;
	  msgpack::adaptor::convert<Ubitrack::Math::Vector<double, 4> >()(o, vec);
      v = Ubitrack::Math::Quaternion::fromVector(vec);
      return o;
  }
};

template<>
struct pack<Ubitrack::Math::Quaternion > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::Quaternion& v) const
  {
      Ubitrack::Math::Vector<double, 4> vec;
      v.toVector(vec);
      o.pack(vec);
      return o;
  }
};

template<>
struct object_with_zone<Ubitrack::Math::Quaternion > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::Quaternion& v) const
  {
      o.type = type::ARRAY;
      o.via.array.size = 4;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      o.via.array.ptr[0] = msgpack::object(v.x(), o.zone);
      o.via.array.ptr[1] = msgpack::object(v.y(), o.zone);
      o.via.array.ptr[2] = msgpack::object(v.z(), o.zone);
      o.via.array.ptr[3] = msgpack::object(v.w(), o.zone);
  }
};



/*
 * Ubitrack::Math::RotationVelocity
 */
#if !defined(MSGPACK_USE_CPP03)
template<>
struct as<Ubitrack::Math::RotationVelocity> {
  Ubitrack::Math::RotationVelocity operator()(msgpack::object const& o) const
  {
      return Ubitrack::Math::RotationVelocity(o.as<Ubitrack::Math::Vector<double, 3> >());
  }
};
#endif

template<>
struct convert<Ubitrack::Math::RotationVelocity > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::RotationVelocity& v) const
  {
	  Ubitrack::Math::Vector<double, 3> vec;
	  msgpack::adaptor::convert<Ubitrack::Math::Vector<double, 3> >()(o, vec);
      v = Ubitrack::Math::RotationVelocity(vec);
      return o;
  }
};

template<>
struct pack<Ubitrack::Math::RotationVelocity > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::RotationVelocity& v) const
  {
      Ubitrack::Math::Vector<double, 3> vec = (Ubitrack::Math::Vector<double, 3>)v;
      o.pack(vec);
      return o;
  }
};

template<>
struct object_with_zone<Ubitrack::Math::RotationVelocity > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::RotationVelocity& v) const
  {
      o.type = type::ARRAY;
      o.via.array.size = 3;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      o.via.array.ptr[0] = msgpack::object(v(0), o.zone);
      o.via.array.ptr[1] = msgpack::object(v(1), o.zone);
      o.via.array.ptr[2] = msgpack::object(v(2), o.zone);
  }
};



/*
 * Ubitrack::Math::Pose
 */
#if !defined(MSGPACK_USE_CPP03)
template<>
struct as<Ubitrack::Math::Pose> {
  Ubitrack::Math::Pose operator()(msgpack::object const& o) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 2) throw msgpack::type_error();
      return Ubitrack::Math::Pose(
              o.via.array.ptr[0].as<Ubitrack::Math::Quaternion>(),
              o.via.array.ptr[1].as<Ubitrack::Math::Vector<double, 3> >()
      );
  }
};
#endif

template<>
struct convert<Ubitrack::Math::Pose > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::Pose& v) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 2) throw msgpack::type_error();
	  Ubitrack::Math::Quaternion quat;
	  Ubitrack::Math::Vector<double, 3> vec;
	  msgpack::adaptor::convert<Ubitrack::Math::Quaternion>()(o.via.array.ptr[0], quat);
	  msgpack::adaptor::convert<Ubitrack::Math::Vector<double, 3> >()(o.via.array.ptr[1], vec);
      v = Ubitrack::Math::Pose(quat, vec);
      return o;
  }
};

template<>
struct pack<Ubitrack::Math::Pose > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::Pose& v) const
  {
      o.pack_array(2);
      o.pack(v.rotation());
      o.pack(v.translation());
      return o;
  }
};

template<>
struct object_with_zone<Ubitrack::Math::Pose > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::Pose& v) const
  {
      o.type = type::ARRAY;
      o.via.array.size = 2;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      o.via.array.ptr[0] = msgpack::object(v.rotation(), o.zone);
      o.via.array.ptr[1] = msgpack::object(v.translation(), o.zone);
  }
};



/*
 * Ubitrack::Math::ErrorVector<T,N>
 */
#if !defined(MSGPACK_USE_CPP03)
template<typename T, std::size_t N>
struct as<Ubitrack::Math::ErrorVector<T, N> > {
  Ubitrack::Math::ErrorVector<T, N> operator()(msgpack::object const& o) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 2) throw msgpack::type_error();
      return Ubitrack::Math::ErrorVector<T, N>(
              o.via.array.ptr[0].as<Ubitrack::Math::Vector<T, N> >(),
              o.via.array.ptr[1].as<Ubitrack::Math::Matrix<T, N, N> >()
      );
  }
};
#endif

template<typename T, std::size_t N>
struct convert<Ubitrack::Math::ErrorVector<T, N> > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::ErrorVector<T, N>& v) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 2) throw msgpack::type_error();
	  Ubitrack::Math::Vector<T, N> vec;
	  Ubitrack::Math::Matrix<T, N, N> cov;
	  msgpack::adaptor::convert<Ubitrack::Math::Vector<T, N> >()(o.via.array.ptr[0], vec);
	  msgpack::adaptor::convert<Ubitrack::Math::Matrix<T, N, N> >()(o.via.array.ptr[1], cov);
      v = Ubitrack::Math::ErrorVector<T, N>(vec, cov);
      return o;
  }
};

template<typename T, std::size_t N>
struct pack<Ubitrack::Math::ErrorVector<T, N> > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::ErrorVector<T, N>& v) const
  {
      o.pack_array(2);
      o.pack(v.value);
      o.pack(v.covariance);
      return o;
  }
};

template<typename T, std::size_t N>
struct object_with_zone<Ubitrack::Math::ErrorVector<T, N> > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::ErrorVector<T, N>& v) const
  {
      o.type = type::ARRAY;
      o.via.array.size = 2;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      o.via.array.ptr[0] = msgpack::object(v.value, o.zone);
      o.via.array.ptr[1] = msgpack::object(v.covariance, o.zone);
  }
};



/*
 * Ubitrack::Math::ErrorPose
 */
#if !defined(MSGPACK_USE_CPP03)
template<>
struct as<Ubitrack::Math::ErrorPose> {
  Ubitrack::Math::ErrorPose operator()(msgpack::object const& o) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 3) throw msgpack::type_error();
      return Ubitrack::Math::ErrorPose(
              o.via.array.ptr[0].as<Ubitrack::Math::Quaternion>(),
              o.via.array.ptr[1].as<Ubitrack::Math::Vector<double, 3> >(),
              o.via.array.ptr[2].as< Ubitrack::Math::Matrix< double, 6, 6 > >()
      );
  }
};
#endif

template<>
struct convert<Ubitrack::Math::ErrorPose > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::ErrorPose& v) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 3) throw msgpack::type_error();
	  Ubitrack::Math::Quaternion quat;
	  Ubitrack::Math::Vector<double, 3> vec;
	  Ubitrack::Math::Matrix<double, 6, 6> cov;
	  msgpack::adaptor::convert<Ubitrack::Math::Quaternion>()(o.via.array.ptr[0], quat);
	  msgpack::adaptor::convert<Ubitrack::Math::Vector<double, 3> >()(o.via.array.ptr[1], vec);
	  msgpack::adaptor::convert<Ubitrack::Math::Matrix<double, 6, 6> >()(o.via.array.ptr[2], cov);
      v = Ubitrack::Math::ErrorPose(quat, vec, cov);
      return o;
  }
};

template<>
struct pack<Ubitrack::Math::ErrorPose > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::ErrorPose& v) const
  {
      o.pack_array(3);
      o.pack(v.rotation());
      o.pack(v.translation());
      o.pack(v.covariance());
      return o;
  }
};

template<>
struct object_with_zone<Ubitrack::Math::ErrorPose > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::ErrorPose& v) const
  {
      o.type = type::ARRAY;
      o.via.array.size = 3;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      o.via.array.ptr[0] = msgpack::object(v.rotation(), o.zone);
      o.via.array.ptr[1] = msgpack::object(v.translation(), o.zone);
      o.via.array.ptr[2] = msgpack::object(v.covariance(), o.zone);
  }
};


/*
 * Ubitrack::Math::CameraIntrinsics<T>
 */
#if !defined(MSGPACK_USE_CPP03)
template<typename T>
struct as<Ubitrack::Math::CameraIntrinsics<T> > {
  Ubitrack::Math::CameraIntrinsics<T> operator()(msgpack::object const& o) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 6) throw msgpack::type_error();
      Ubitrack::Math::CameraIntrinsics<T> result;
      result.calib_type = static_cast<typename Ubitrack::Math::CameraIntrinsics<T>::CalibType>(o.via.array.ptr[0].as<int>());
      result.dimension = o.via.array.ptr[1].as<Ubitrack::Math::Vector< std::size_t, 2 > >();
      result.matrix = o.via.array.ptr[2].as<Ubitrack::Math::Matrix< T, 3, 3 > >();
      result.radial_size = o.via.array.ptr[3].as<std::size_t>();
      result.radial_params = o.via.array.ptr[4].as<Ubitrack::Math::Vector< T, 6 > >();
      result.tangential_params = o.via.array.ptr[5].as<Ubitrack::Math::Vector< T, 2 > >();
      return result;
  }
};
#endif

template<typename T>
struct convert<Ubitrack::Math::CameraIntrinsics<T> > {
  msgpack::object const& operator()(msgpack::object const& o, Ubitrack::Math::CameraIntrinsics<T>& v) const
  {
      if (o.type!=msgpack::type::ARRAY) throw msgpack::type_error();
      if (o.via.array.size != 6) throw msgpack::type_error();
	  int calib_type;
	  msgpack::adaptor::convert<int>()(o.via.array.ptr[0], calib_type);
      v.calib_type = static_cast<typename Ubitrack::Math::CameraIntrinsics<T>::CalibType>(calib_type);
	  msgpack::adaptor::convert<Ubitrack::Math::Vector<std::size_t, 2> >()(o.via.array.ptr[1], v.dimension);
	  msgpack::adaptor::convert<Ubitrack::Math::Matrix< T, 3, 3 > >()(o.via.array.ptr[2], v.matrix);
	  msgpack::adaptor::convert<std::size_t >()(o.via.array.ptr[3], v.radial_size);
	  msgpack::adaptor::convert<Ubitrack::Math::Vector< T, 6 > >()(o.via.array.ptr[4], v.radial_params);
	  msgpack::adaptor::convert<Ubitrack::Math::Vector< T, 2 > >()(o.via.array.ptr[5], v.tangential_params);
      return o;
  }
};

template<typename T>
struct pack<Ubitrack::Math::CameraIntrinsics<T> > {
  template<typename Stream>
  msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, const Ubitrack::Math::CameraIntrinsics<T>& v) const
  {
      o.pack_array(6);
      o.pack((int)v.calib_type);
      o.pack(v.dimension);
      o.pack(v.matrix);
      o.pack(v.radial_size);
      o.pack(v.radial_params);
      o.pack(v.tangential_params);
      return o;
  }
};

template<typename T>
struct object_with_zone<Ubitrack::Math::CameraIntrinsics<T> > {
  void operator()(msgpack::object::with_zone& o, const Ubitrack::Math::CameraIntrinsics<T>& v) const
  {
      o.type = type::ARRAY;
      o.via.array.size = 6;
      o.via.array.ptr = static_cast<msgpack::object*>(
              o.zone.allocate_align(sizeof(msgpack::object)*o.via.array.size));
      o.via.array.ptr[0] = msgpack::object((int)v.calib_type, o.zone);
      o.via.array.ptr[1] = msgpack::object(v.dimension, o.zone);
      o.via.array.ptr[2] = msgpack::object(v.matrix, o.zone);
      o.via.array.ptr[3] = msgpack::object(v.radial_size, o.zone);
      o.via.array.ptr[4] = msgpack::object(v.radial_params, o.zone);
      o.via.array.ptr[5] = msgpack::object(v.tangential_params, o.zone);
  }
};



} // namespace adaptor
} // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack


#endif //UBITRACK_MSGPACKSERIALIZER_H

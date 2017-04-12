
#include <utSerialization/MsgpackSerializer.h>
#include <utMeasurement/Measurement.h>
#include <utUtil/Exception.h>

#include <string>
#include <fstream>
#include <vector>
#include <algorithm> //std::transform

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack;
using namespace Ubitrack::Serialization;

#ifdef HAVE_MSGPACK

template< typename T >
void testSerializeSimple(const T& data)
{
	// Serialize
    msgpack::sbuffer ostream;
	MsgpackArchive::serialize(ostream, data);

	// Deserialize
	T result;

    MsgpackArchive::deserialize(ostream, result);

	BOOST_CHECK_EQUAL(data, result);
}

template< typename T >
void testSerializeSimpleVector(const std::vector<T>& data)
{
    // Serialize
    msgpack::sbuffer ostream;
    MsgpackArchive::serialize(ostream, data);

    // Deserialize
    std::vector<T> result;

    MsgpackArchive::deserialize(ostream, result);

    BOOST_CHECK_EQUAL(data.size(), result.size());
    for (std::size_t i=0; i<data.size(); ++i) {
        BOOST_CHECK_EQUAL(data.at(i), result.at(i));

    }
}


template< typename T >
void testSerializeMeasurement(const Measurement::Measurement< T >& data)
{
    // Serialize
    msgpack::sbuffer ostream;
    MsgpackArchive::serialize(ostream, data);

    // Deserialize
    Measurement::Measurement< T > result = Measurement::Measurement< T >( 0, boost::shared_ptr< T >( new T() ) );

    MsgpackArchive::deserialize(ostream, result);

    BOOST_CHECK_EQUAL(data.time(), result.time());
    BOOST_CHECK_EQUAL(*data, *result);
}

template< typename T >
void testSerializeMeasurementVector(const Measurement::Measurement< std::vector< T > >& data)
{
    // Serialize
    msgpack::sbuffer ostream;
    MsgpackArchive::serialize(ostream, data);

    // Deserialize
    Measurement::Measurement< std::vector< T > > result = Measurement::Measurement< std::vector< T > >( 0, boost::shared_ptr< std::vector< T > >( new std::vector<T>() ) );

    MsgpackArchive::deserialize(ostream, result);

    BOOST_CHECK_EQUAL(data.time(), result.time());
    BOOST_CHECK_EQUAL(data->size(), result->size());
    for (std::size_t i=0; i<data->size(); ++i) {
        BOOST_CHECK_EQUAL(data->at(i), result->at(i));
    }

}

template<typename T>
std::vector<T> make_vector_simple(const T& v, int count) {
    std::vector<T> result(count);
    for (std::size_t i=0; i<count; ++i) {
        result.push_back(v);
    }
    return result;
}

template<typename T>
Ubitrack::Measurement::Measurement< std::vector< T > > make_vector_measurement(unsigned long long ts, const T& v, int count) {
    std::vector<T> result(count);
    for (std::size_t i=0; i<count; ++i) {
        result.push_back(v);
    }
    return Ubitrack::Measurement::Measurement< std::vector< T > >(ts, result);
}



void testSerializeMultiple()
{

    // Math::Scalar<int>
    Math::Scalar<int> v_scalari(22);
    // Math::Scalar<double>
    Math::Scalar<double> v_scalard(22.33);
    // Vector
    Math::Vector<double, 3> v_vec3 = randomVector<double, 3>(5.0);
    // Quaternion
    Math::Quaternion v_quat = randomQuaternion();
    // Pose
    Math::Pose v_pose(randomQuaternion(), randomVector<double, 3>(5.0));
    // Matrix3x3
    Math::Matrix<double, 3, 3> v_mat33;
    randomMatrix(v_mat33);
    // Matrix4x4
    Math::Matrix<double, 4, 4> v_mat44;
    randomMatrix(v_mat44);

    msgpack::sbuffer buffer;
    msgpack::packer<msgpack::sbuffer> pk(&buffer);

    // serialize all elements
    MsgpackArchive::serialize(pk, v_scalari);
    MsgpackArchive::serialize(pk, v_scalard);
    MsgpackArchive::serialize(pk, v_vec3);
    MsgpackArchive::serialize(pk, v_quat);
    MsgpackArchive::serialize(pk, v_pose);
    MsgpackArchive::serialize(pk, v_mat33);
    MsgpackArchive::serialize(pk, v_mat44);

    // transfer to destination
    msgpack::unpacker pac;
    pac.reserve_buffer(buffer.size());
    // simulate transfer
    memcpy(pac.buffer(), buffer.data(), buffer.size() );

    pac.buffer_consumed(buffer.size());

    // Deserialize elements
    Math::Scalar<int> r_scalari;
    Math::Scalar<double> r_scalard;
    Math::Vector<double, 3> r_vec3;
    Math::Quaternion r_quat;
    Math::Pose r_pose;
    Math::Matrix<double, 3, 3> r_mat33;
    Math::Matrix<double, 4, 4> r_mat44;

    MsgpackArchive::deserialize(pac, r_scalari);
    BOOST_CHECK_EQUAL(v_scalari, r_scalari);
    MsgpackArchive::deserialize(pac, r_scalard);
    BOOST_CHECK_EQUAL(v_scalard, r_scalard);
    MsgpackArchive::deserialize(pac, r_vec3);
    BOOST_CHECK_EQUAL(v_vec3, r_vec3);
    MsgpackArchive::deserialize(pac, r_quat);
    BOOST_CHECK_EQUAL(v_quat, r_quat);
    MsgpackArchive::deserialize(pac, r_pose);
    BOOST_CHECK_EQUAL(v_pose, r_pose);
    MsgpackArchive::deserialize(pac, r_mat33);
    BOOST_CHECK_EQUAL(v_mat33, r_mat33);
    MsgpackArchive::deserialize(pac, r_mat44);
    BOOST_CHECK_EQUAL(v_mat44, r_mat44);
}



#endif // HAVE_MSGPACK


void TestMsgpack()
{

#ifdef HAVE_MSGPACK
    // Test simple data types

    // Math::Scalar<int>
    Math::Scalar<int> v_scalari(22);
    testSerializeSimple(v_scalari);

    std::vector< Math::Scalar<int> > vec_scalari = make_vector_simple(v_scalari, 5);
    testSerializeSimpleVector(vec_scalari);

    // Math::Scalar<double>
    Math::Scalar<double> v_scalard(22.33);
    testSerializeSimple(v_scalard);

    std::vector< Math::Scalar<double> > vec_scalard = make_vector_simple(v_scalard, 5);
    testSerializeSimpleVector(vec_scalard);

    // Vector
    Math::Vector<double, 3> v_vec3 = randomVector<double, 3>(5.0);
    testSerializeSimple(v_vec3);

    std::vector< Math::Vector<double, 3> > vec_vec3 = make_vector_simple(v_vec3, 5);
    testSerializeSimpleVector(vec_vec3);

    // Quaternion
    Math::Quaternion v_quat = randomQuaternion();
    testSerializeSimple(v_quat);

    std::vector< Math::Quaternion > vec_quat = make_vector_simple(v_quat, 5);
    testSerializeSimpleVector(vec_quat);

    // Pose
    Math::Pose v_pose(randomQuaternion(), randomVector<double, 3>(5.0));
    testSerializeSimple(v_pose);

    std::vector< Math::Pose > vec_pose = make_vector_simple(v_pose, 5);
    testSerializeSimpleVector(vec_pose);

    // Matrix3x3
    Math::Matrix<double, 3, 3> v_mat33;
    randomMatrix(v_mat33);
    testSerializeSimple(v_mat33);

    std::vector< Math::Matrix<double, 3, 3> > vec_mat33 = make_vector_simple(v_mat33, 5);
    testSerializeSimpleVector(vec_mat33);

    // Matrix4x4
    Math::Matrix<double, 4, 4> v_mat44;
    randomMatrix(v_mat44);
    testSerializeSimple(v_mat44);

    std::vector< Math::Matrix<double, 4, 4> > vec_mat44 = make_vector_simple(v_mat44, 5);
    testSerializeSimpleVector(vec_mat44);




    // Test Measurements
	Measurement::Timestamp ts = Measurement::now();

    // Button
    Measurement::Button m_button(ts, v_scalari);
    testSerializeMeasurement(m_button);

    Measurement::ButtonList vec_button = make_vector_measurement(ts, v_scalari, 5);
    testSerializeMeasurementVector(vec_button);

    // Distance
	Measurement::Distance m_distance(ts, v_scalard);
	testSerializeMeasurement(m_distance);

    Measurement::DistanceList vec_distance = make_vector_measurement(ts, v_scalard, 5);
    testSerializeMeasurementVector(vec_distance);

    // Position
    Measurement::Position m_pos(ts, v_vec3);
    testSerializeMeasurement(m_pos);

    Measurement::PositionList vec_pos = make_vector_measurement(ts, v_vec3, 5);
    testSerializeMeasurementVector(vec_pos);

    // Rotation
    Measurement::Rotation m_quat(ts, v_quat);
    testSerializeMeasurement(m_quat);

    // Pose
    Measurement::Pose m_pose(ts, v_pose);
    testSerializeMeasurement(m_pose);

    Measurement::PoseList vec_posem = make_vector_measurement(ts, v_pose, 5);
    testSerializeMeasurementVector(vec_posem);

    // Matrix3x3
    Measurement::Matrix3x3 m_mat33(ts, v_mat33);
    testSerializeMeasurement(m_mat33);

    // Matrix4x4
    Measurement::Matrix4x4 m_mat44(ts, v_mat44);
    testSerializeMeasurement(m_mat44);


    // test serializing multiple objects in astream
    testSerializeMultiple();
    
#endif // HAVE_MSGPACK
}



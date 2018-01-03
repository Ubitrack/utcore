
#include <utSerialization/BoostArchiveSerializer.h>
#include <utMeasurement/Measurement.h>
#include <utUtil/Exception.h>

#include <string>
#include <fstream>
#include <algorithm> //std::transform

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack;
using namespace Ubitrack::Serialization;

template< typename T >
void testSerializeSimpleText(const T& data)
{
	// Serialize
	std::ostringstream ostream;
	boost::archive::text_oarchive out_archive( ostream );
	BoostArchive::serialize(out_archive, data);
    // terminate string
    out_archive << "\n";
	std::string serialized_data(ostream.str());

	// Deserialize
	T result;

	std::istringstream istream(serialized_data);
	boost::archive::text_iarchive in_archive( istream );
	BoostArchive::deserialize(in_archive, result);

	BOOST_CHECK_EQUAL(data, result);
}

template< typename T >
void testSerializeSimpleBinary(const T& data)
{
	// Serialize
	std::ostringstream ostream;
	boost::archive::binary_oarchive out_archive( ostream );
	BoostArchive::serialize(out_archive, data);

	std::string serialized_data(ostream.str());

	// Deserialize
	T result;

	typedef boost::iostreams::basic_array_source<char> Device;
	boost::iostreams::stream_buffer<Device> istream((char*)serialized_data.data(), serialized_data.size());
	boost::archive::binary_iarchive in_archive( istream );
    BoostArchive::deserialize(in_archive, result);

	BOOST_CHECK_EQUAL(data, result);
}


template< typename T >
void testSerializeMeasurementText(const Measurement::Measurement< T >& data)
{
    // Serialize
    std::ostringstream ostream;
    boost::archive::text_oarchive out_archive( ostream );
    BoostArchive::serialize(out_archive, data);
    // terminate string
    out_archive << "\n";
    std::string serialized_data(ostream.str());

    // Deserialize
    Measurement::Measurement< T > result = Measurement::Measurement< T >( 0, boost::shared_ptr< T >( new T() ) );

    std::istringstream istream(serialized_data);
    boost::archive::text_iarchive in_archive( istream );
    BoostArchive::deserialize(in_archive, result);

    BOOST_CHECK_EQUAL(data.time(), result.time());
    BOOST_CHECK_EQUAL(*data, *result);
}

template< typename T >
void testSerializeMeasurementBinary(const Measurement::Measurement< T >& data)
{
    // Serialize
    std::ostringstream ostream;
    boost::archive::binary_oarchive out_archive( ostream );
    BoostArchive::serialize(out_archive, data);

    std::string serialized_data(ostream.str());

    // Deserialize
    Measurement::Measurement< T > result = Measurement::Measurement< T >( 0, boost::shared_ptr< T >( new T() ) );

    typedef boost::iostreams::basic_array_source<char> Device;
    boost::iostreams::stream_buffer<Device> istream((char*)serialized_data.data(), serialized_data.size());
    boost::archive::binary_iarchive in_archive( istream );
    BoostArchive::deserialize(in_archive, result);

    BOOST_CHECK_EQUAL(data.time(), result.time());
    BOOST_CHECK_EQUAL(*data, *result);
}


void testSerializeBinaryMultiple()
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


    // Serialize
    std::ostringstream ostream;
    boost::archive::binary_oarchive out_archive( ostream );
    // serialize all elements
    BoostArchive::serialize(out_archive, v_scalari);
    BoostArchive::serialize(out_archive, v_scalard);
    BoostArchive::serialize(out_archive, v_vec3);
    BoostArchive::serialize(out_archive, v_quat);
    BoostArchive::serialize(out_archive, v_pose);
    BoostArchive::serialize(out_archive, v_mat33);
    BoostArchive::serialize(out_archive, v_mat44);

    std::string serialized_data(ostream.str());

    // Deserialize elements
    Math::Scalar<int> r_scalari;
    Math::Scalar<double> r_scalard;
    Math::Vector<double, 3> r_vec3;
    Math::Quaternion r_quat;
    Math::Pose r_pose;
    Math::Matrix<double, 3, 3> r_mat33;
    Math::Matrix<double, 4, 4> r_mat44;

    typedef boost::iostreams::basic_array_source<char> Device;
    boost::iostreams::stream_buffer<Device> istream((char*)serialized_data.data(), serialized_data.size());
    boost::archive::binary_iarchive in_archive( istream );

    BoostArchive::deserialize(in_archive, r_scalari);
    BOOST_CHECK_EQUAL(v_scalari, r_scalari);
    BoostArchive::deserialize(in_archive, r_scalard);
    BOOST_CHECK_EQUAL(v_scalard, r_scalard);
    BoostArchive::deserialize(in_archive, r_vec3);
    BOOST_CHECK_EQUAL(v_vec3, r_vec3);
    BoostArchive::deserialize(in_archive, r_quat);
    BOOST_CHECK_EQUAL(v_quat, r_quat);
    BoostArchive::deserialize(in_archive, r_pose);
    BOOST_CHECK_EQUAL(v_pose, r_pose);
    BoostArchive::deserialize(in_archive, r_mat33);
    BOOST_CHECK_EQUAL(v_mat33, r_mat33);
    BoostArchive::deserialize(in_archive, r_mat44);
    BOOST_CHECK_EQUAL(v_mat44, r_mat44);
}



void TestBoostArchive()
{
    // Test simple data types

    // Math::Scalar<int>
    Math::Scalar<int> v_scalari(22);
    testSerializeSimpleText(v_scalari);
    testSerializeSimpleBinary(v_scalari);

    // Math::Scalar<double>
    Math::Scalar<double> v_scalard(22.33);
    testSerializeSimpleText(v_scalard);
    testSerializeSimpleBinary(v_scalard);

    // Vector
    Math::Vector<double, 3> v_vec3 = randomVector<double, 3>(5.0);
    testSerializeSimpleText(v_vec3);
    testSerializeSimpleBinary(v_vec3);

    // Quaternion
    Math::Quaternion v_quat = randomQuaternion();
    testSerializeSimpleText(v_quat);
    testSerializeSimpleBinary(v_quat);

    // Pose
    Math::Pose v_pose(randomQuaternion(), randomVector<double, 3>(5.0));
    testSerializeSimpleText(v_pose);
    testSerializeSimpleBinary(v_pose);

    // Matrix3x3
    Math::Matrix<double, 3, 3> v_mat33;
    randomMatrix(v_mat33);
    testSerializeSimpleText(v_mat33);
    testSerializeSimpleBinary(v_mat33);

    // Matrix4x4
    Math::Matrix<double, 4, 4> v_mat44;
    randomMatrix(v_mat44);
    testSerializeSimpleText(v_mat44);
    testSerializeSimpleBinary(v_mat44);



    // Test Measurements
    Measurement::Timestamp ts = Measurement::now();

    // Button
    Measurement::Button m_button(ts, v_scalari);
    testSerializeMeasurementText(m_button);
    testSerializeMeasurementBinary(m_button);

    // Distance
    Measurement::Distance m_distance(ts, v_scalard);
    testSerializeMeasurementText(m_distance);
    testSerializeMeasurementBinary(m_distance);

    // Position
    Measurement::Position m_pos(ts, v_vec3);
    testSerializeMeasurementText(m_pos);
    testSerializeMeasurementBinary(m_pos);

    // Rotation
    Measurement::Rotation m_quat(ts, v_quat);
    testSerializeMeasurementText(m_quat);
    testSerializeMeasurementBinary(m_quat);

    // Pose
    Measurement::Pose m_pose(ts, v_pose);
    testSerializeMeasurementText(m_pose);
    testSerializeMeasurementBinary(m_pose);

    // Matrix3x3
    Measurement::Matrix3x3 m_mat33(ts, v_mat33);
    testSerializeMeasurementText(m_mat33);
    testSerializeMeasurementBinary(m_mat33);

    // Matrix4x4
    Measurement::Matrix4x4 m_mat44(ts, v_mat44);
    testSerializeMeasurementText(m_mat44);
    testSerializeMeasurementBinary(m_mat44);

    // test serializing multiple objects in astream
    testSerializeBinaryMultiple();
}



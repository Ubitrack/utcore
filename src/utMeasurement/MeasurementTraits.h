//
// Created by Ulrich Eck on 2019-08-05.
//

#ifndef UBITRACK_CORE_MEASUREMENTTRAITS_H
#define UBITRACK_CORE_MEASUREMENTTRAITS_H

#include <string>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include "utMeasurement/Measurement.h"

namespace Ubitrack {
    namespace Measurement {
        namespace Traits {


// Enumeration of all available element data types
            enum class DataType {
                DataTypeUndefined = 0,
                DataTypeInteger,
                DataTypeUnsignedLong,
                DataTypeDouble
            };


// Enumeration of all available Types for introspection
            enum class MeasurementType {
                MeasurementTypeUndefined = 0,
                MeasurementTypeScalar,
                MeasurementTypeVector2,
                MeasurementTypeVector3,
                MeasurementTypeVector4,
                MeasurementTypeVector8,
                MeasurementTypeQuaternion,
                MeasurementTypeMatrix3x3,
                MeasurementTypeMatrix3x4,
                MeasurementTypeMatrix4x4,
                MeasurementTypePose,
                MeasurementTypeErrorPose,
                MeasurementTypeErrorVector2,
                MeasurementTypeErrorVector3,
                MeasurementTypeRotationVelocity,
                MeasurementTypeCameraIntrinsics,
                MeasurementTypeImage // Forward declaration as we cannot extend an enumeration later
            };


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
 * \brief A vector datatype is one whose size is not constant, i.e. it has variable-length arrays or strings
 */
            template<typename M>
            struct IsVector : public FalseType {
            };

            // there should be a generic way to do this for all Math::TYPES that are wrapped in a vector ..
            // what's wrong with this??
//            template<>
//            struct IsVector<Measurement< std::vector < M::value_type > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::Scalar< int > > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::Scalar< double > > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::Scalar< unsigned long > > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::Pose > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::Vector< double, 2 > > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::Vector< double, 3 > > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::ErrorPose > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::ErrorVector< double, 2 > > > > : public TrueType {};

            template<>
            struct IsVector<Measurement< std::vector < Math::ErrorVector< double, 3 > > > > : public TrueType {};



            template< typename T >
            struct MeasurementTypeTraits
            {
                bool isFixedType() const
                { return IsVector<typename boost::remove_reference<typename boost::remove_const<T>::type>::type>::value; }

                DataType getDataType() const
                { return DataType ::DataTypeUndefined; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeUndefined; }

            };


            // button
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Button >
             {
                DataType getDataType() const
                { return DataType ::DataTypeInteger; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeScalar; }
            };

            // distance
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Distance >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeScalar; }
            };

            // Vector2
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Position2D >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeVector2; }
            };

            // Vector3
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Position >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeVector3; }
            };

            // Vector4
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Vector4D >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeVector4; }
            };

            // Quaternion
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Rotation >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeQuaternion; }
            };

            // Matrix3x3
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Matrix3x3 >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeMatrix3x3; }
            };

            // Matrix3x4
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Matrix3x4 >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeMatrix3x4; }
            };

            // Matrix4x4
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Matrix4x4 >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeMatrix4x4; }
            };

            // Pose
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::Pose >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypePose; }
            };

            // ErrorPose
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::ErrorPose >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeErrorPose; }
            };

            // ErrorPosition2
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::ErrorPosition2 >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeErrorVector2; }
            };

            // ErrorPosition3
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::ErrorPosition >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeErrorVector3; }
            };

            // RotationVelocity
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::RotationVelocity >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeRotationVelocity; }
            };

            // CameraIntrinsics
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::CameraIntrinsics >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeCameraIntrinsics; }
            };


            // multiple measurements

            // Vector of Buttons
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::ButtonList >
            {
                DataType getDataType() const
                { return DataType ::DataTypeInteger; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeScalar; }
            };

            // Vector of Distances
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::DistanceList >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeScalar; }
            };

            // Vector of IDs
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::IDList >
            {
                DataType getDataType() const
                { return DataType ::DataTypeUnsignedLong; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeScalar; }
            };

            // Vector of Poses
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::PoseList >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypePose; }
            };

            // Vector of Position2
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::PositionList2 >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeVector2; }
            };

            // Vector of Position
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::PositionList >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeVector3; }
            };

            // Vector of ErrorPoses
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::ErrorPoseList >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeErrorPose; }
            };

            // Vector of ErrorPositionList2
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::ErrorPositionList2 >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeErrorVector2; }
            };

            // Vector of ErrorPosition
            template<>
            struct MeasurementTypeTraits< Ubitrack::Measurement::ErrorPositionList >
            {
                DataType getDataType() const
                { return DataType ::DataTypeDouble; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::MeasurementTypeErrorVector3; }
            };

        } // Traits
    } // Measurement
} // Ubitrack


#endif //UBITRACK_CORE_MEASUREMENTTRAITS_H

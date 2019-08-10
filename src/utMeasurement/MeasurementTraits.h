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


// Enumeration of all available Types for introspection
            enum class MeasurementType {
                Undefined = 0,
                ScalarInt,
                ScalarDouble,
                ScalarUnsignedLong,
                Vector2,
                Vector3,
                Vector4,
                Vector8,
                Quaternion,
                Matrix3x3,
                Matrix3x4,
                Matrix4x4,
                Pose,
                ErrorPose,
                ErrorVector2,
                ErrorVector3,
                RotationVelocity,
                CameraIntrinsics,
                Image // Forward declaration as we cannot extend an enumeration later
            };

            std::ostream& operator<<( std::ostream& s, const MeasurementType& m )
            {
                switch (m) {
                    case Undefined:
                        s << "Undefined";
                        break;
                    case ScalarInt:
                        s << "ScalarInt";
                        break;
                    case ScalarDouble:
                        s << "ScalarDouble";
                        break;
                    case ScalarUnsignedLong:
                        s << "ScalarUnsignedLong";
                        break;
                    case Vector2:
                        s << "Vector2";
                        break;
                    case Vector3:
                        s << "Vector3";
                        break;
                    case Vector4:
                        s << "Vector4";
                        break;
                    case Vector8:
                        s << "Vector8";
                        break;
                    case Quaternion:
                        s << "Quaternion";
                        break;
                    case Matrix3x3:
                        s << "Matrix3x3";
                        break;
                    case Matrix3x4:
                        s << "Matrix3x4";
                        break;
                    case Matrix4x4:
                        s << "Matrix4x4";
                        break;
                    case Pose:
                        s << "Pose";
                        break;
                    case ErrorPose:
                        s << "ErrorPose";
                        break;
                    case ErrorVector2:
                        s << "ErrorVector2";
                        break;
                    case ErrorVector3:
                        s << "ErrorVector3";
                        break;
                    case RotationVelocity:
                        s << "RotationVelocity";
                        break;
                    case CameraIntrinsics:
                        s << "CameraIntrinsics";
                        break;
                    case Image:
                        s << "Image";
                        break;
                }
                return s;
            }


            template< typename T >
            struct MeasurementTypeToEnumTraits
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Undefined; }
            };


            template< MeasurementType MT, bool ISV = false >
            struct MeasurementEnumToTypeTraits
            {
            };

            // Measurement to Enum Traits

            // button
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Button >
             {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ScalarInt; }
            };

            // distance
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Distance >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ScalarDouble; }
            };

            // Vector2
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Position2D >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Vector2; }
            };

            // Vector3
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Position >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Vector3; }
            };

            // Vector4
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Vector4D >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Vector4; }
            };

            // Vector8
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Vector8D >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Vector8; }
            };

            // Quaternion
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Rotation >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Quaternion; }
            };

            // Matrix3x3
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Matrix3x3 >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Matrix3x3; }
            };

            // Matrix3x4
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Matrix3x4 >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Matrix3x4; }
            };

            // Matrix4x4
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Matrix4x4 >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Matrix4x4; }
            };

            // Pose
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::Pose >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Pose; }
            };

            // ErrorPose
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::ErrorPose >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ErrorPose; }
            };

            // ErrorPosition2
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::ErrorPosition2 >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ErrorVector2; }
            };

            // ErrorPosition3
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::ErrorPosition >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ErrorVector3; }
            };

            // RotationVelocity
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::RotationVelocity >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::RotationVelocity; }
            };

            // CameraIntrinsics
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::CameraIntrinsics >
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::CameraIntrinsics; }
            };


            // multiple measurements

            // Vector of Buttons
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::ButtonList >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ScalarInt; }
            };

            // Vector of Distances
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::DistanceList >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ScalarDouble; }
            };

            // Vector of IDs
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::IDList >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ScalarUnsignedLong; }
            };

            // Vector of Poses
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::PoseList >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Pose; }
            };

            // Vector of Position2
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::PositionList2 >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Vector2; }
            };

            // Vector of Position
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::PositionList >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Vector3; }
            };

            // Vector of ErrorPoses
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::ErrorPoseList >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ErrorPose; }
            };

            // Vector of ErrorPositionList2
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::ErrorPositionList2 >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ErrorVector2; }
            };

            // Vector of ErrorPosition
            template<>
            struct MeasurementTypeToEnumTraits< Ubitrack::Measurement::ErrorPositionList >
            {
                bool isFixedType() const
                { return false; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::ErrorVector3; }
            };



            // Enum to Measurement Traits

            // Button
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ScalarInt, false>
            {
                typedef Ubitrack::Measurement::Button value_type;
            };

            // Distance
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ScalarDouble, false>
            {
                typedef Ubitrack::Measurement::Distance value_type;
            };

            // Position2D
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Vector2, false>
            {
                typedef Ubitrack::Measurement::Position2D value_type;
            };

            // Position
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Vector3, false>
            {
                typedef Ubitrack::Measurement::Position value_type;
            };

            // Vector4d
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Vector4, false>
            {
                typedef Ubitrack::Measurement::Vector4D value_type;
            };

            // Vector8d
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Vector8, false>
            {
                typedef Ubitrack::Measurement::Vector8D value_type;
            };

            // Matrix3x3
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Matrix3x3, false>
            {
                typedef Ubitrack::Measurement::Matrix3x3 value_type;
            };

            // Matrix3x4
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Matrix3x4, false>
            {
                typedef Ubitrack::Measurement::Matrix3x4 value_type;
            };

            // Matrix4x4
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Matrix4x4, false>
            {
                typedef Ubitrack::Measurement::Matrix4x4 value_type;
            };

            // Rotation
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Quaternion, false>
            {
                typedef Ubitrack::Measurement::Rotation value_type;
            };

            // Pose
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Pose, false>
            {
                typedef Ubitrack::Measurement::Pose value_type;
            };

            // ErrorPosition2
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ErrorVector2, false>
            {
                typedef Ubitrack::Measurement::ErrorPosition2 value_type;
            };

            // ErrorPosition
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ErrorVector3, false>
            {
                typedef Ubitrack::Measurement::ErrorPosition value_type;
            };

            // ErrorPose
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ErrorPose, false>
            {
                typedef Ubitrack::Measurement::ErrorPose value_type;
            };

            // RotationVelocity
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::RotationVelocity, false>
            {
                typedef Ubitrack::Measurement::RotationVelocity value_type;
            };

            // CameraIntrinsics
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::CameraIntrinsics, false>
            {
                typedef Ubitrack::Measurement::CameraIntrinsics value_type;
            };

            // Vector of Buttons
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ScalarInt, true>
            {
                typedef Ubitrack::Measurement::ButtonList value_type;
            };

            // Vector of Distances
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ScalarDouble, true>
            {
                typedef Ubitrack::Measurement::DistanceList value_type;
            };

            // Vector of IDs
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ScalarUnsignedLong, true>
            {
                typedef Ubitrack::Measurement::IDList value_type;
            };

            // Vector of Poses
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Pose, true>
            {
                typedef Ubitrack::Measurement::PoseList value_type;
            };

            // Vector of Position2
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Vector2, true>
            {
                typedef Ubitrack::Measurement::PositionList2 value_type;
            };

            // Vector of Position
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::Vector3, true>
            {
                typedef Ubitrack::Measurement::PositionList value_type;
            };

            // Vector of ErrorPoses
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ErrorPose, true>
            {
                typedef Ubitrack::Measurement::ErrorPoseList value_type;
            };

            // Vector of ErrorPositionList2
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ErrorVector2, true>
            {
                typedef Ubitrack::Measurement::ErrorPositionList2 value_type;
            };

            // Vector of ErrorPosition
            template<>
            struct MeasurementEnumToTypeTraits<MeasurementType::ErrorVector3, true>
            {
                typedef Ubitrack::Measurement::ErrorPositionList value_type;
            };

        } // Traits
    } // Measurement
} // Ubitrack


#endif //UBITRACK_CORE_MEASUREMENTTRAITS_H

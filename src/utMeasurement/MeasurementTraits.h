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



/**
 * \brief Base type for compile-type true/false tests.  Compatible with Boost.MPL.  classes inheriting from this type
 * are \b false values.
 */
            struct FalseType {
                static const bool value = false;
                typedef FalseType type;
            };

            template< typename T >
            struct MeasurementTypeToEnumTraits
            {
                bool isFixedType() const
                { return true; }

                MeasurementType getMeasurementType() const
                { return MeasurementType ::Undefined; }
            };



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

        } // Traits
    } // Measurement
} // Ubitrack


#endif //UBITRACK_CORE_MEASUREMENTTRAITS_H

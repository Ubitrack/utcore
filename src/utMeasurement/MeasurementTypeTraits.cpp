//
// Created by Ulrich Eck on 2019-08-10.
//

#include "MeasurementTraits.h"

using namespace Ubitrack::Measurement::Traits;

std::ostream& operator<<( std::ostream& s, const MeasurementType& m )
{
    switch (m) {
        case MeasurementType::Undefined:
            s << "Undefined";
            break;
        case MeasurementType::ScalarInt:
            s << "ScalarInt";
            break;
        case MeasurementType::ScalarDouble:
            s << "ScalarDouble";
            break;
        case MeasurementType::ScalarUnsignedLong:
            s << "ScalarUnsignedLong";
            break;
        case MeasurementType::Vector2:
            s << "Vector2";
            break;
        case MeasurementType::Vector3:
            s << "Vector3";
            break;
        case MeasurementType::Vector4:
            s << "Vector4";
            break;
        case MeasurementType::Vector8:
            s << "Vector8";
            break;
        case MeasurementType::Quaternion:
            s << "Quaternion";
            break;
        case MeasurementType::Matrix3x3:
            s << "Matrix3x3";
            break;
        case MeasurementType::Matrix3x4:
            s << "Matrix3x4";
            break;
        case MeasurementType::Matrix4x4:
            s << "Matrix4x4";
            break;
        case MeasurementType::Pose:
            s << "Pose";
            break;
        case MeasurementType::ErrorPose:
            s << "ErrorPose";
            break;
        case MeasurementType::ErrorVector2:
            s << "ErrorVector2";
            break;
        case MeasurementType::ErrorVector3:
            s << "ErrorVector3";
            break;
        case MeasurementType::RotationVelocity:
            s << "RotationVelocity";
            break;
        case MeasurementType::CameraIntrinsics:
            s << "CameraIntrinsics";
            break;
        case MeasurementType::Image:
            s << "Image";
            break;
    }
    return s;
}

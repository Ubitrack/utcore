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
 * @ingroup tracking_algorithms
 * @file
 * Functions for Absolute Orientation (3D-3D-Pose-Estimation)
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_ABSOLUTE_ORIENTATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_ABSOLUTE_ORIENTATION_H_INCLUDED__



#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>

#include <vector>

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation3D3D {

/**
 * @ingroup tracking_algorithms
 * @internal
 * Calculate the absolute orientation problem.
 * This method calculates the pose between two coorinate frames
 * as specified by corresponding pairs of 3D points. (Absolute
 * Orientation Problem).
 * The implementation is according to Horn "Closed-form solution
 * of the absolute orientation using unit quaternions"
 * (J. Optical Soc. of America A, Vol. 4, page 629, 1987),
 * except the quaternion is determined by the eigenvalue problem and
 * not by Ferrari's method.
 *
 * @param left vector of 3D points in the left coordinate frame.
 * @param right vector of 3D points in the right coordinate frame.
 * @return pose that describes the transformation of the left coordinate frame into the right coordinate frame.
 * @throws Util::Exception if lapack is not available, different number of samples for left and right side are given or the matrix N only has non-positive eigenvalues.
 */
UBITRACK_EXPORT Math::Pose calculateAbsoluteOrientation ( const std::vector< Math::Vector< double, 3 > >& left,
														  const std::vector< Math::Vector< double, 3 > >& right);

/**
 * @ingroup tracking_algorithms
 * @brief Calculates a solution to the 3D-3D pose estiataion problem, also
 * know as the Absolute Orientation problem.
 *
 * This algrotihm calculates the pose between two coordinate frames
 * as specified by corresponding pairs of 3D points. (Absolute
 * Orientation Problem). The implementation is based on "Closed-form
 * solution of the absolute orientation using unit quaternions"
 * from Horn ( @cite horn1987closed ).
 * @verbatim
@article{horn1987closed,
  title={Closed-form solution of absolute orientation using unit quaternions},
  author={Horn, Berthold KP},
  journal={JOSA A},
  volume={4},
  number={4},
  pages={629--642},
  year={1987},
  publisher={Optical Society of America}
} @endverbatim
 *
 * A major difference to the original publication is the solution
 * using eigenvalues instead of Ferrari's method.
 *
 * Example use case:\n
 * @code
 * std::vector< Vector3d > points3d1; // <- should be filled with 3D object points
 * std::vector< Vector3d > points3d2; // <- should be filled with corresponding 3D points in another basis
 * Math::Pose pose; // <- will be filled with the transformation
 * estimatePose6D_3D3D( points3d1, pose, points3d1 );
 * @endcode
 *
 * @attention : There is a version of this function ovlerloaded with \c float instead of \c double parameters.
 *
 * @param points3dA vector of \b 3D \b points in the left coordinate frame.
 * @return pose the \b pose describes the transformation of the right coordinate frame into the left coordinate frame.
 * @param points3dB vector of \b 3D \b points in the right coordinate frame.
 * @return a flag that signs if the algorithm has succesfully determined a solution.
 */
UBITRACK_EXPORT bool estimatePose6D_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Pose& pose
										, const std::vector< Math::Vector3d >& points3dB );

/** 
 * @internal
 * @brief overloaded function \c estimatePose6D_3D3D with \c float parameters.
 *
 * For further information on this algorithm see estimatePose6D_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Pose& pose, const std::vector< Math::Vector3d >& points3dB );
 */
UBITRACK_EXPORT bool estimatePose6D_3D3D( const std::vector< Math::Vector3f >& points3dA, Math::Pose& pose
										, const std::vector< Math::Vector3f >& points3dB );

/** 
 * @brief This algorithm estimates the rotation between two coordinate frames.
 *
 * This algorithm is a solution to the subproblem of the Aboslute Orientation problem.
 *
 * For further information on this algorithm see estimatePose6D_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Pose& pose, const std::vector< Math::Vector3d >& points3dB );
 *
 * @attention : There is a version of this function ovlerloaded with \c float instead of \c double parameters and some others with \e quaternion as return value.
 *
 * @param points3dA vector of \b 3D \b points in the left coordinate frame.
 * @return mat describes the orientation from the right coordinate frame into the left coordinate frame as a \e 3-by-3 \e matrix.
 * @param points3dB vector of \b 3D \b points in the right coordinate frame.
 * @return a flag that signs if the algorithm has succesfully determined a solution.
 */
UBITRACK_EXPORT bool estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Matrix3x3d& mat
										, const std::vector< Math::Vector3d >& points3dB );

/** 
 * @internal
 * @brief overloaded function \c estimateRotation_3D3D with \c float parameters.
 *
 * For further information on this algorithm see estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Matrix3x3d& mat, const std::vector< Math::Vector3d >& points3dB );
 */
UBITRACK_EXPORT bool estimateRotation_3D3D( const std::vector< Math::Vector3f >& points3dA, Math::Matrix3x3f& mat
										, const std::vector< Math::Vector3f >& points3dB );
										

/** 
 * @brief overloaded function \c estimateRotation_3D3D with \c double parameters and a \e quaternion instead of a \e matrix as a return value.
 *
 * For further information on this algorithm see estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Matrix3x3d& mat, const std::vector< Math::Vector3d >& points3dB );
 */
UBITRACK_EXPORT bool estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Quaternion& quat
										, const std::vector< Math::Vector3d >& points3dB );
										
/** 
 * @internal
 * @brief overloaded function \c estimateRotation_3D3D with \c float parameters and a \e quaternion instead of a \e matrix as a return value.
 *
 * For further information on this algorithm see estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Matrix3x3d& mat, const std::vector< Math::Vector3d >& points3dB );
 */
UBITRACK_EXPORT bool estimateRotation_3D3D( const std::vector< Math::Vector3f >& points3dA, Math::Quaternion& quat
										, const std::vector< Math::Vector3f >& points3dB );

/** 
 * @brief This algorithm estimates the scale between two points clouds.
 *
 *
 * For further information on this algorithm see estimatePose6D_3D3D( const std::vector< Math::Vector3d >& points3dA, Math::Pose& pose, const std::vector< Math::Vector3d >& points3dB );
 *
 * @attention : There is a version of this function ovlerloaded with \c float instead of \c double parameters.
 *
 * @param m_left vector of \b 3D \b points in the left coordinate frame.
 * @param m_right vector of coressponding \b 3D \b points in another coordinate frame with a differen scale.
 * @return the scale that can be applied to scale the second points cloud to similar lengths as the first one. 
 */
UBITRACK_EXPORT double estimateScale_3D3D ( const std::vector< Math::Vector3d >& points3dA
									, const std::vector< Math::Vector3d >& points3dB );

/** 
 * @internal
 * @brief overloaded function \c estimateScale_3D3D with \c float parameters instead of \c double .
 *
 * For further information on this algorithm see estimateScale_3D3D( const std::vector< Math::Vector3d >& points3dA, const std::vector< Math::Vector3d >& points3dB );
 */
UBITRACK_EXPORT float estimateScale_3D3D ( const std::vector< Math::Vector3f >& points3dA
									, const std::vector< Math::Vector3f >& points3dB );
	
/**
 * @internal function object version of calculateAbsoluteOrientation for RANSAC etc.
 */
template< class T >
class EstimateAbsoluteOrientation
{
public:
	void operator()( Math::Pose& result, const std::vector< Math::Vector< T, 3 > >& pointsA, 
		const std::vector< Math::Vector< T, 3 > >& pointsB ) const
	{
		Math::Pose pose( calculateAbsoluteOrientation( pointsA, pointsB ) );
		result = pose;
	}
};

/**
 * @internal function object to evalute an absolute orientation for RANSAC etc.
 */
template< class T >
class EvaluateAbsoluteOrientation
{
public:
	/**
	 * @internal computes euclidean distance of transformed point to original point
	 */
	T operator()( const Math::Pose &p, const Math::Vector< T, 3 > &a, const Math::Vector< T, 3 > &b ) const
	{
		Math::Vector< T, 3 > c = p * a;
		return boost::numeric::ublas::norm_2( b - c );
	}
};

// the following should go into Stochastic folder namespace:
// /** 
 // * @brief Returns the Root Mean Square of two point sets that have been transformed by a pose
 // * @param pose The transformation from the left to the right point cloud
 // * @param left The left point cloud
 // * @param right The right point cloud
 // * @returns Rms of corresponding point distances
 // */
// template< class T >
// T computeRms( const Math::Pose& pose, const std::vector< Math::Vector< T, 3 > >& left,
    // const std::vector< Math::Vector< T, 3 > >& right )
// {
    // if ( left.size() != right.size() )
        // UBITRACK_THROW( "Invalid input list sizes for computeRms" );
    // double rms = 0.0;
    // for ( size_t i = 0; i < left.size(); ++i ) {
        // rms +=  boost::numeric::ublas::norm_2( pose * left.at( i ) - right.at( i ) );
    // }
    // return rms;
// }

} } } // namespace Ubitrack::Algorithm::PoseEstimation3D3D

#endif // HAVE_LAPACK

#endif

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
#include <utMath/Scalar.h>
#include <vector>

namespace Ubitrack { namespace Calibration {

UBITRACK_EXPORT Math::Scalar< double > calculateAbsoluteOrientationScale ( const std::vector< Math::Vector< 3, double > >& m_left,
														         const std::vector< Math::Vector< 3, double > >& m_right);

/**
 * @ingroup tracking_algorithms
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
UBITRACK_EXPORT Math::Pose calculateAbsoluteOrientation ( const std::vector< Math::Vector< 3, double > >& left,
														  const std::vector< Math::Vector< 3, double > >& right);

/**
 * @ingroup tracking_algorithms
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
 * @param leftBegin pointer at the beginning of the 3D vectors in the left coordinate frame.
 * @param leftEnd pointer at one after the last of the 3D vectors in the left coordinate frame.
 * @param rightBegin pointer at the beginning of the 3D vectors in the right coordinate frame.
 * @param rightEnd pointer at one after the last of the 3D vectors in the left coordinate frame.
 * @return pose that describes the transformation of the left coordinate frame into the right coordinate frame.
 * @throws Util::Exception if lapack is not available, different number of samples for left and right side are given or the matrix N only has non-positive eigenvalues.
 */
UBITRACK_EXPORT Math::Pose calculateAbsoluteOrientation ( const Math::Vector< 3, double >* leftBegin, const Math::Vector< 3, double >* leftEnd,
														  const Math::Vector< 3, double >* rightBegin, const Math::Vector< 3, double >* rightEnd );

/**
 * @ingroup tracking_algorithms
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
 * @param leftBegin std::vector iterator that points at the beginning of the 3D vectors in the left coordinate frame.
 * @param leftEnd std::vector iterator that points at one after the last of the 3D vectors in the left coordinate frame.
 * @param rightBegin std::vector iterator that points at the beginning of the 3D vectors in the right coordinate frame.
 * @param rightEnd std::vector iterator that points at one after the last of the 3D vectors in the left coordinate frame.
 * @return pose that describes the transformation of the left coordinate frame into the right coordinate frame.
 * @throws Util::Exception if lapack is not available, different number of samples for left and right side are given or the matrix N only has non-positive eigenvalues.
 */
UBITRACK_EXPORT Math::Pose calculateAbsoluteOrientation ( const std::vector< Math::Vector< 3, double > >::iterator& leftBegin,
														  const std::vector< Math::Vector< 3, double > >::iterator& leftEnd,
														  const std::vector< Math::Vector< 3, double > >::iterator& rightBegin,
														  const std::vector< Math::Vector< 3, double > >::iterator& rightEnd );

/**
 * function object version of calculateAbsoluteOrientation for RANSAC etc.
 */

template< class T >
class EstimateAbsoluteOrientation
{
public:
	void operator()( Math::Pose& result, const std::vector< Math::Vector< 3, T > >& pointsA, 
		const std::vector< Math::Vector< 3, T > >& pointsB ) const
	{
		Math::Pose pose( calculateAbsoluteOrientation( pointsA, pointsB ) );
		result = pose;
	}
};

/**
 * function object to evalute an absolute orientation for RANSAC etc.
 */
template< class T >
class EvaluateAbsoluteOrientation
{
public:
	/**
	 * computes euclidean distance of transformed point to original point
	 */
	T operator()( const Math::Pose &p, const Math::Vector< 3, T > &a, const Math::Vector< 3, T > &b ) const
	{
		Math::Vector< 3, T > c = p * a;
		return boost::numeric::ublas::norm_2( b - c );
	}
};

/** 
 * @brief Returns the Root Mean Square of two point sets that have been transformed by a pose
 * @param pose The transformation from the left to the right point cloud
 * @param left The left point cloud
 * @param right The right point cloud
 * @returns Rms of corresponding point distances
 */
template< class T >
T computeRms( const Math::Pose& pose, const std::vector< Math::Vector< 3, T > >& left,
    const std::vector< Math::Vector< 3, T > >& right )
{
    if ( left.size() != right.size() )
        UBITRACK_THROW( "Invalid input list sizes for computeRms" );
    double rms = 0.0;
    for ( size_t i = 0; i < left.size(); ++i ) {
        rms +=  boost::numeric::ublas::norm_2( pose * left.at( i ) - right.at( i ) );
    }
    return rms;
}

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

#endif

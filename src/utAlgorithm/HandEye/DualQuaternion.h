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
 * A Dual Quaternion solution to the hand-eye calibration problem.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGORITHM_HANDEYE_DUAL_QUATERNION_CALIBRATION_H_INCLUDED__
#define __UBITRACK_ALGORITHM_HANDEYE_DUAL_QUATERNION_CALIBRATION_H_INCLUDED__


#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/Pose.h>

#include <vector>

namespace Ubitrack { namespace Algorithm { namespace HandEye {

/**
 * @ingroup calibration tracking_algorithms
 * @brief An algorithm to determine a solution to the classic \b Hand-Eye
 * \b Calibration problem, based on given \b 6D \b pose correspondences.
 
 * This algorithm estimates a \b pose from given \b 6D \b pose correspondences.
 * This problem is well known from robotics research but is also of special 
 * interest for Augmented Reality scenarios. Among all the many solutions
 * that can be found to this problem the implementation of this solution is 
 * based on the the publication 'Hand-Eye Calibration Using Dual Quaternions'
 * from Konstantinos Daniilidis in 1999 ( @cite daniilidis1999hand ).
 *
 * The hand eye calibration can be seen as a solution to determine the a priori
 * unkown pose \b p that specifies a spatial transformation from one coordinate
 * frame to another one which are rigidly connected to each other. Several
 * observations (at least three) in each coordinate frame are necessary
 * to determine a solution to the hand-eye problem.\n
 * If @f$ a_i * p * b_i \f$ describes this spatial transformation, typically \n
 * @f$ a_i @f$ are n poses in a camera coordinate frame, specifying the pose from
 * the camera (eye) to an observed target, that is usually rigidly connected to
 * a robot or similar. \n
 * @f$ b_i @f$ are n poses in a robots' coordinate frame, specifying the pose from
 * the robots' base to the robots' hand (holding the camera). \n
 * The number of n(n-1)/2 distinct pose pairs are used to determine the
 * solution using their pose differences.
 *
 * @verbatim
 @article{daniilidis1999hand,
  title={Hand-eye calibration using dual quaternions},
  author={Daniilidis, Konstantinos},
  journal={The International Journal of Robotics Research},
  volume={18},
  number={3},
  pages={286--298},
  year={1999},
  publisher={SAGE Publications}
} @endverbatim
 *
 * Example use case:\n
 @code
 std::vector< Pose > posesA; // <- should be filled with poses in one coordinate system (eye)
 std::vector< Pose > posesB; // <- should be filled with corresponding poses from another coordinate system (hand)
 Math::Pose pose; // <- will be filled with a solution, if there is one.
 estimatePose6D_6D6D( posesA, pose, posesB );
 @endcode
 *
 * @attention : Other versions might occur in future, this algorithm is still under development.
 * @note : Also have a look at \b select_6DPoses and \b generate_relative_6DPoses wihich might be of interest in future versions
 *
 * @param eyes \b 6D \b poses in the \b 1st coordinate system.
 * @param pose the pose, if a solution can be found.
 * @param hands corresponding \b 6D \b poses in the \b 2nd coordinate system.
 * @return a flag that signs if the algorithm has succesfully determined a solution.
 */
UBITRACK_EXPORT bool estimatePose6D_6D6D(  const std::vector< Math::Pose >& eyes, Math::Pose& pose,
	const std::vector< Math::Pose >& hands );

}}} // namespace Ubitrack::Algorithm::HandEye

#endif // HAVE_LAPACK

#endif //__UBITRACK_ALGORITHM_HANDEYE_DUAL_QUATERNION_CALIBRATION_H_INCLUDED__

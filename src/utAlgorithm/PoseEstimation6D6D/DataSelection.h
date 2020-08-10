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

#ifndef __UBITRACK_ALGORITHM_HANDEYE_DATA_SELECTION_H_INCLUDED__
#define __UBITRACK_ALGORITHM_HANDEYE_DATA_SELECTION_H_INCLUDED__


#include <utCore.h>
#include <utMath/Pose.h> //includes quaternion
#include <utMath/Vector.h>
#include <utMath/Util/RotationCast.h>
#include <utMath/Geometry/container_traits.h>
#include <utMath/Stochastic/identity_iterator.h>
#include <utMath/Stochastic/k_means.h> // copy_probability, k_means

#include <vector>

namespace Ubitrack{ namespace Math {


/// @todo maybe rename to geodesic distance later?
template< typename pose_type >
struct RotationDistance
{
};

template< typename T >
struct RotationDistance< Math::Vector< T, 6 > >
{
	typedef Math::Vector< T, 6 > pose_type;
	T operator() ( const Math::Vector< T, 6 >& pose1, const Math::Vector< T, 6 >& pose2 ) const
	{
		const T x = pose1[ 0 ] - pose2[ 0 ];
		const T y = pose1[ 1 ] - pose2[ 1 ];
		const T z = pose1[ 2 ] - pose2[ 2 ];
		return std::sqrt( x*x+y*y+z*z );
	}
};

/**
 * @internal Functor that casts a pose of one type into a pose of another one type.
 *
 * This functor can be applied whenever pose data available in
 * one representation needs to be transferred into another representation.
 * Of course the information is not changed, despite numerical errors, and
 * this functor can be integrated into templated functions to limit 
 * the rewriting of algorithm code.
 *
 * @tparam pose_type defines the type of the pose to cast to.
 */
template< typename pose_type >
struct pose_cast
{
	typedef pose_type return_type;
	
	template< typename pose_in >
	return_type operator()( const pose_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_MATCHING_POSE_CONVERSION );
		return pose_type();
	}
};

/** @internal specialization of functor to cast in to the \b ubitrack default pose \b format */
template< >
struct pose_cast< Math::Pose  >
{
	typedef Math::Pose result_type;

	template< typename pose_in >
	result_type operator()( const pose_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_MATCHING_POSE_CONVERSION );
		return result_type();
	}

	template< typename T >
	result_type operator()( const Math::Vector< T, 6 > &pose ) const
	{
		const Math::Vector< T, 3 > rotationAxis( pose[ 0 ], pose[ 1 ], pose[ 2 ] );
		const Math::Quaternion quatRot( Math::Util::RotationCast< Math::Quaternion >()( rotationAxis ) );
		const Math::Vector< double, 3 > translation( pose[ 3 ], pose[ 4 ], pose[ 5 ] );

		return Math::Pose( quatRot, translation );
	}
	
	template< typename T >
	result_type operator()( const Math::Vector< T, 8 > &dualQuat ) const
	{
		const Math::Quaternion quatRot( dualQuat[ 1 ], dualQuat[ 2 ], dualQuat[ 3 ], dualQuat[ 0 ] );
		const Math::Quaternion qprime( dualQuat[ 5 ], dualQuat[ 6 ], dualQuat[ 7 ], dualQuat[ 4 ] );
		const Math::Quaternion q_conj( -dualQuat[ 1 ], -dualQuat[ 2 ], -dualQuat[ 3 ], dualQuat[ 0 ] );
		const Math::Quaternion t_final = qprime*q_conj;

		return Math::Pose( quatRot, Math::Vector< double, 3 >( 2*t_final.x(), 2*t_final.y(), 2*t_final.z() ) );
	}
};

/**
 * @internal functor specialization to cast into a 7-vector pose
 * representation that consists of a 4-vector angle/axis rotation
 * and a 3-vector translation part.
 * briefly: angle, rotAxisX, rotAxisY, rotAxisZ, transX, transY, transZ
 */
template< typename T >
struct pose_cast< Math::Vector< T, 7 >  >
{
	typedef Math::Vector< T, 7 > result_type;

	template< typename rotation_in >
	result_type operator()( const rotation_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_MATCHING_POSE_CONVERSION );
		return result_type();
	}


	result_type operator()( const Math::Pose &pose ) const
	{
		Math::Vector< T, 7 > result;
		
		Math::Vector< T, 4 > rot = Math::Util::RotationCast< Math::Vector< T, 4 > >()( pose.rotation() );
		result[ 0 ] = rot[ 0 ];
		result[ 1 ] = rot[ 1 ];
		result[ 2 ] = rot[ 2 ];
		result[ 3 ] = rot[ 3 ];
		result[ 4 ] = pose.translation()[ 0 ];
		result[ 5 ] = pose.translation()[ 1 ];
		result[ 6 ] = pose.translation()[ 2 ];
		return result;
	}
};

/**
 * @internal functor specialization to cast into a 6-vector pose
 * representation that consists of a 3-vector rotation axis
 * where the angle is encoded as the length of the rotation axis
 * and a 3-vector translation part.
 * briefly: rotAxisX, rotAxisY, rotAxisZ, transX, transY, transZ
 */
template< typename T >
struct pose_cast< Math::Vector< T, 6 >  >
{
	typedef Math::Vector< T, 6 > result_type;

	template< typename pose_in >
	result_type operator()( const pose_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_MATCHING_POSE_CONVERSION );
		return result_type();
	}


	result_type operator()( const Math::Pose &pose ) const
	{
		Math::Vector< T, 6 > result;
		
		Math::Vector< T, 3 > rot = Math::Util::RotationCast< Math::Vector< T, 3 > >()( pose.rotation() );
		result[ 0 ] = rot[ 0 ];
		result[ 1 ] = rot[ 1 ];
		result[ 2 ] = rot[ 2 ];
		result[ 3 ] = pose.translation()[ 0 ];
		result[ 4 ] = pose.translation()[ 1 ];
		result[ 5 ] = pose.translation()[ 2 ];
		return result;
	}
	
	template< typename Type >
	result_type operator()( const Math::Vector< Type, 8 >& dualQuat ) const
	{
		const Math::Pose pose = pose_cast< Math::Pose >()( dualQuat );
		const Math::Vector< T, 3 > rotAxis = Math::Util::RotationCast< Math::Vector< T, 3 > >()( pose.rotation() );
		Math::Vector< T, 6 > result;
		result[ 0 ] = rotAxis[ 0 ];
		result[ 1 ] = rotAxis[ 1 ];
		result[ 2 ] = rotAxis[ 2 ];
		result[ 3 ] = pose.translation()[ 0 ];
		result[ 4 ] = pose.translation()[ 1 ];
		result[ 5 ] = pose.translation()[ 2 ];
		return result;
	}
};

/**
 * @internal functor specialization to cast into a 8-vector pose
 * representation that consists of a 4-vector quaternion encoding
 * the rotation and a 4-vector quaternion encoding the translation.
 *
 * briefly: rotQuatW, rotQuatX, rotQuatY, rotQuatZ
 * , transQuatW, transQuatX, transQuatY, transQuatZ,
 */
template< typename T >
struct pose_cast< Math::Vector< T, 8 > >
{
	typedef typename Math::Vector< T, 8 > result_type;
	
	Math::Vector< T, 8 > operator()( const Math::Pose& pose ) const
	{
		const T qw = pose.rotation().w();
		const T qx = pose.rotation().x();
		const T qy = pose.rotation().y();
		const T qz = pose.rotation().z();
		
		const T tx = pose.translation()[ 0 ];
		const T ty = pose.translation()[ 1 ];
		const T tz = pose.translation()[ 2 ];
		
		// dual-quaternion return value: (q | q')
		// order: q  :=  (qw  | qx  qy  qz )
		// order: q' :=  (q'w | q'x q'y q'z)
		Math::Vector< T, 8 > dualQuat;
		
		// first, the easy part: :)
		// quaternion goes into quaternion part (q) 
		dualQuat( 0 ) = qw;
		dualQuat( 1 ) = qx;
		dualQuat( 2 ) = qy;
		dualQuat( 3 ) = qz;
				
		// second, the tricky part:
		// translation t goes into quaternion dual part (q' == aka q prime )
		// how does it work?
		// you take translation t and make a quaternion q_t out of it,
		// assuming zero for the real part an (t/2) as the imaginary parts.
		// you apply quaternion multiplication to q_t with q from the right side
		// and done. 
		// summary: q' := q( t/2, 0) * q
		//
		// technically it is done in two parts
		// part 1: inner product of imaginary part from q and t
		// -(1/2) * (q_xyz' *t) 
		dualQuat( 4 ) = -0.5* ( tx*qx + ty*qy + tz*qz );
		
		// part 2:
		// (1/2) * (cross_product(t, q_xyz) + q_w*t)
		dualQuat( 5 ) = 0.5* ( ( ty*qz - tz*qy ) + qw*tx );
		dualQuat( 6 ) = 0.5* ( ( tz*qx - tx*qz ) + qw*ty );
		dualQuat( 7 ) = 0.5* ( ( tx*qy - ty*qx ) + qw*tz );
		
		return dualQuat;
	}
};

/**
 * @internal functor that transforms rotations that
 * are of arbitrary alignment into the same hemisphere.
 *
 * The functor accepts any type of pose or rotation data.
 * naja, it should accept, actually it does not. :)
 */
template< typename pose_type, bool positive >
struct hemisphere_alignment
{
	typedef pose_type result_type;

	template< typename pose_in >
	result_type operator()( const pose_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_MATCHING_SPECIALIZATION_AVALIABLE );
		return result_type();
	}
};

template< typename T >
struct hemisphere_alignment< Math::Vector< T, 6 >, true > 
{
	typedef Math::Vector< T, 6 > result_type;

	template< typename pose_in >
	result_type operator()( const pose_in & pose )const
	{
		UBITRACK_STATIC_ASSERT( (false), NO_MATCHING_POSE_AVAILABLE );
		return result_type();
	}
	
	Math::Vector< T, 6 > operator()( const Math::Vector< T, 6 > & pose ) const
	{
		Math::Vector< T, 6 > result( pose );
		const T x = pose[ 0 ];
		const T y = pose[ 1 ];
		const T z = pose[ 2 ];
		const T angle = std::sqrt( x*x+y*y+z*z );
		
		if( angle > 3.1415926535897932384626433832795 )
		{
			const T new_angle = (2*3.1415926535897932384626433832795) - angle;
			const T ratio = new_angle/angle;
			result[ 0 ] *= ratio;
			result[ 1 ] *= ratio;
			result[ 2 ] *= ratio;
		}
		return result;
	}
};

/** 
 * @internal Functor that calculates a relative pose from two given poses.
 *
 * The relative pose as the result of this functor is calculated from two
 * given absolute poses. The relative pose can be uses as an input parameter
 * for presumably all types of hand-eye calibration algorithms. 
 *
 * In the example of the dual quaternion approach to the hand-eye calibration
 * the resulting pose would be labled as @f$ a @f$ .
 *
 * @tparam pose_type defining the type of the resulting output pose
 * ( e.g. Pose, DualQuaternion, etc.)
 * @tparam inv signs if the forward ( @f$ a @f$ )or the backward
 * ( @f$ b @f$ )difference pose is calculated.
 */
 
template< typename pose_type, bool inv >
struct relative_pose
{
	typedef pose_type return_type;
	
	return_type operator()( const Math::Pose& pose1, const Math::Pose& pose2 ) const
	{
		const Math::Pose pose = (~pose2) * pose1;
		return pose_cast< return_type >()( pose );
	}
};

/** @internal this functor specialization calculates a relative pose for the backward constellation, see aboove */
template< typename pose_type >
struct relative_pose< pose_type, false >
{
	typedef pose_type return_type;

	return_type operator()( const Math::Pose& pose1, const Math::Pose& pose2 ) const
	{
		const Math::Pose pose = (pose2) * (~pose1);
		return pose_cast< return_type >()( pose );
	}
};

namespace Util {
/** 
 * @internal
 * @brief 
 * This function computes a result among the input values from first to last using the given binary operator
 *
 * This template function is inspired from http://www.cplusplus.com/reference/numeric/adjacent_difference/
 */
template< class InputIterator, class OutputIterator, class BinaryOperation >
OutputIterator adjacent_difference ( InputIterator first, InputIterator last, OutputIterator result, BinaryOperation binary_op )
{
	if (first!=last)
	{
		typename std::iterator_traits< InputIterator >::value_type val,prev;
		prev = *first;
		while(++first!=last)
		{
			val = *first;
			*++result++ = binary_op( val,prev );
			prev = val;
		}
		++result;
	}
	return result;
}

} // namespace Ubitrack::Math::Util

}} //namesapace Ubitrack::Math;


namespace Ubitrack { namespace Algorithm { namespace PoseEstimation6D6D {

/**
 * @brief this structure supports the process of selecting best matches of relative pose measurements for a hand eye calibration.
 *
 */
template< typename T >	
// struct UBITRACK_EXPORT DataSelection
struct DataSelection
{
	/// type of the underlying precision (e.g. \c double or \c float )
	typedef T precision_type;
	
	/// at the moment the selection process works only with axis angle
	typedef Math::Vector< precision_type, 6 > pose_type;
	
	/// sequence container type used for storing the prepared relative pose measurements
	typedef std::vector< pose_type > RelPoseListType;
	
	/// type of index sequence that 
	typedef std::vector< std::size_t > IndexListType;
	
	
	RelPoseListType comparisonPoseList;
	
	RelPoseListType clusterPoseList;
	
	IndexListType indexList;

	/// this function alings all orientation of pose measurements in the upper hemisphere
	template< typename PoseType >
	pose_type operator( )( const PoseType& pose_in ) const
	{
		const pose_type pose = Math::pose_cast< pose_type >()( pose_in );
		return Math::hemisphere_alignment< pose_type, true >()( pose );
	}
	
	template< typename PoseType, template< typename Type, typename = std::allocator< Type > > class container_type >
	inline void resetComparisonPoses( const std::size_t n, const container_type< PoseType >& relPosesIn )
	{
		// const std::size_t n = std::distance( relPosesIn.begin(), relPosesIn.end() );
		updateComparisonPoses( relPosesIn.begin(), relPosesIn.end() , comparisonPoseList );
		
		// update the cluster of best matches of relative pose measurements
		clusterPoseList.clear();
		clusterPoseList.reserve( n );
		resetCluster( n, comparisonPoseList.begin(), comparisonPoseList.end(), std::back_inserter( clusterPoseList ) );
		// update indices of relative poses that are nearest to the cluster centers
		resetIndices( comparisonPoseList.begin(), comparisonPoseList.end(), clusterPoseList.begin(), clusterPoseList.end() );
	}
	
	/// function returns the best matches of a relative pose measurements depending on indices
	template< typename PoseType, template< typename Type, typename = std::allocator< Type > > class container_type, typename OutputIterator >
	void getSelection( const container_type< PoseType >& relPosesIn, OutputIterator itSelectedPoseOut )
	{
		std::vector< size_t >::const_iterator itIndex = indexList.begin();
		const std::vector< size_t >::const_iterator itIndexEnd = indexList.end();
		
		for( ; itIndex != itIndexEnd; ++itIndex)
		{
			typename container_type< PoseType >::const_iterator it = relPosesIn.begin();
			std::advance( it, *itIndex );
			*itSelectedPoseOut++ = *it;
		}
	}
	
protected:

	template< typename InputIterator >
	inline void updateComparisonPoses( const InputIterator itBegin, const InputIterator itEnd , RelPoseListType& relPoseList ) const
	{
		const std::size_t nIn = std::distance ( itBegin, itEnd );
		const std::size_t nIs = std::distance ( relPoseList.begin(), relPoseList.end() );
		assert( nIn <= nIs );
		
		InputIterator it = itBegin;
		std::advance( it, nIs );
		std::transform( it, itEnd, std::back_inserter( relPoseList ), DataSelection< T >() );
	}
	
	template< typename PoseType, template< typename Type, typename = std::allocator< Type > > class container_type >
	inline void updateComparisonPoses( const container_type< PoseType >& relPosesIn )
	{
		updateComparisonPoses( relPosesIn.begin(), relPosesIn.end() , comparisonPoseList );
	}
	

	template< typename InputIterator, typename OutputIterator >
	inline void resetCluster( const std::size_t n_cluster, const InputIterator itBeginValues, const InputIterator itEndValues, OutputIterator itOut )
	{
		const std::size_t newCluster = Math::Stochastic::copy_probability( itBeginValues, itEndValues, n_cluster, itOut, Math::RotationDistance< pose_type >() );
	}
	
	/// function that determines the indices of the relative pose measurements depending on a cluster indicating best matches and a sequence of special aligned poses.
	template< typename InputIterator1, typename InputIterator2 >
	void resetIndices( const InputIterator1 itBeginValues, const InputIterator1 itEndValues, const InputIterator2 itBeginCluster, const InputIterator2 itEndCluster )
	{
		const std::size_t n_values = std::distance( itBeginValues, itEndValues );
		const std::size_t n_cluster = std::distance( itBeginCluster, itEndCluster );
		
		// get the list of indices signing the membership of a single pose to a corresponding mean
		IndexListType indices;
		indices.reserve( n_values );
		Math::Stochastic::k_means( itBeginValues, itEndValues, itBeginCluster, itEndCluster, std::back_inserter( indices ), Math::RotationDistance< pose_type >() );
		
		// set the list of indices to selected poses
		indexList.assign( n_cluster, 0 ); // indexList.reserve( n_cluster );
		std::vector< precision_type > min_distance( n_cluster, 10000 ); // <- value should be greater than maximal expected distance (==pi)
		
		InputIterator1 itPose = itBeginValues;
		IndexListType::const_iterator itIndex = indices.begin();
		const IndexListType::const_iterator itIndexEnd = indices.end();
		
		for( std::size_t i = 0; itIndex != itIndexEnd; ++itIndex, *itPose , ++i )
		{
			// choose codebook entry to lookup corresponding pose measurement
			InputIterator2 itCodebook = itBeginCluster;	
			std::advance( itCodebook, *itIndex );
			
			const precision_type d = Math::RotationDistance< pose_type >()( *itCodebook, *itPose );
			
			typename std::vector< precision_type >::iterator itDist = min_distance.begin();
			std::advance( itDist, *itIndex );
			if( d < *itDist )
			{
				*itDist = d;
				IndexListType::iterator it = indexList.begin();
				std::advance( it, *itIndex );
				*it = i;
			}
		}
	}
};

template< bool use_all_pairs, bool direction, typename InputIterator, typename OutputIterator >
void generate_relative_pose6D_impl( const InputIterator itBegin, const InputIterator itEnd, OutputIterator itOut )
{
	typedef typename Ubitrack::Util::container_traits< OutputIterator >::value_type desired_pose_type;
	
	const std::size_t n_in = std::distance( itBegin, itEnd );
	assert( n_in > 2 );
	
	//depending on the input generate the corresponding amount of relative movements
	if( !use_all_pairs )
		Math::Util::adjacent_difference ( itBegin, itEnd, itOut , Ubitrack::Math::relative_pose< desired_pose_type, direction >() );
	else
	{
		InputIterator itPose = itBegin;
		for( ;itPose != itEnd; )
		{
			Ubitrack::Util::identity< const Math::Pose > id( *itPose );
			std::advance( itPose, 1u );
			std::transform( itPose , itEnd, id.cbegin(), itOut, Ubitrack::Math::relative_pose< desired_pose_type, direction >() );
		}
	}
}

/**
 * @ingroup calibration tracking_algorithms
 * @brief An algorithm to determine adequate relative poses as an input 
 * to  \b Hand-Eye \b Calibration approaches.
 
 * This algorithm selects \b relative \b pose correspondences from given
 * \b 6D \b relative \b pose correspondences.
 * The quality of an Hand-Eye calibration depends on the input data.
 * Several publications describes the effects of the input data onto the 
 * resulting final pose. As an example relative poses with a small angle
 * should be rejected for the calibration procedure as it leads to numerical
 * instability or singular results.
 * The method implemented here was introduced from Schmidt & Niermann 2008
 * and uses a clustering method for automatic pose selection ( @cite schmidt2008data ).
 * Actually the method for clustering could be exchanged with other methods later.
 *
 * @verbatim
@article{schmidt2008data,
  title={Data selection for hand-eye calibration: a vector quantization approach},
  author={Schmidt, Jochen and Niemann, Heinrich},
  journal={The International Journal of Robotics Research},
  volume={27},
  number={9},
  pages={1027--1053},
  year={2008},
  publisher={SAGE Publications}
} @endverbatim
 *
 * Example use case:\n
 @code
 std::vector< Pose > relPosesA; // <- should be filled with relative poses in one coordinate system (eye)
 std::vector< Pose > relPosesB; // <- should be filled with corresponding relative poses from another coordinate system (hand)
 std::vector< Pose > selectedRelPosesA; //  <- will be filled with n (n = select) relative poses selected from the input data
 std::vector< Pose > selectedRelPosesB; //  <- will be filled with n (n = select) corresponding relative poses from the input data
 
 Math::Pose pose; // <- will be filled with a solution, if there is one.
 estimatePose6D_6D6D( relPosesA, relPosesB, n, selectedRelPosesA, selectedRelPosesB );
 @endcode
 *
 * @attention : This version implements method 3.6 from the publication. Other methods might come later...
 *
 * @param eyes \b relative \b 6D \b poses in the \b 1st coordinate system.
 * @param hands corresponding \b relative \b 6D \b poses in the \b 2nd coordinate system.
 * @param select amount of relative poses to select
 * @param eyesOut contains n ( n=select) \b relative \b 6D \b poses from the from the input data
 * @param handsOut contains n ( n=select) corresponding \b relative \b 6D \b poses from the input data
 */
UBITRACK_EXPORT void select_6DPoses( const std::vector< Math::Pose >& eyes, const std::vector< Math::Pose >& hands
	, const std::size_t select
	, std::vector< Math::Pose >& eyesOut, std::vector< Math::Pose >& handsOut );
	

UBITRACK_EXPORT void generate_relative_6DPoses( const std::vector< Math::Pose >& poses
	, std::vector< Math::Vector< double, 8 > >& relativePoses, bool direction_flag );

}}} // namespace Ubitrack::Algorithm::PoseEstimation6D6D

#endif //__UBITRACK_ALGORITHM_HANDEYE_DATA_SELECTION_H_INCLUDED__

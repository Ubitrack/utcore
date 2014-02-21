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
 * Implementation of data selection methods for the various 
 * hand eye calibration approaches. Contanis also lots of
 * functions for common data manpulation in the calibration
 * procedures.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */


#include "HandEyeDataSelection.h"
#include <utMath/Stochastic/k_means.h> // codebook generation

#include <algorithm> //std::transform


namespace Ubitrack { namespace Calibration {

namespace {

template< bool use_all_pairs, typename T, typename InputIterator, typename OutputIterator >
void pose6D_selection_impl( const InputIterator itBeginEye, const InputIterator itEndEye,
	const InputIterator itBeginHand, const InputIterator itEndHand,
	const std::size_t n_select, OutputIterator itEye, OutputIterator itHand )
{
	typedef Ubitrack::Math::Vector< T, 6 > pose_type;
	typedef typename Ubitrack::Util::container_traits< OutputIterator >::value_type output_type;
	
	const std::size_t n_in = std::distance( itBeginEye, itEndEye );
	const std::size_t n_2 = std::distance( itBeginHand, itEndHand );
	assert( n_in == n_2 );
	assert( n_in > 2 );
	
	const std::size_t n = n_in-1;
	const std::size_t m = use_all_pairs ? (n*n_in)/2 : n;
	
	std::vector< pose_type > dualA; // <- paper from Daniilidis uses a and b to specify dual quaternions
	{ // set the first dual quaternion from input differences
		dualA.reserve( m );
		generate_relative_pose6D_impl< use_all_pairs, true >( itBeginEye, itEndEye, std::back_inserter( dualA ) );
	}
	
	std::vector< pose_type > dualB;  // <- paper from Daniilidis uses a and b to specify dual quaternions
	{
		dualB.reserve( m );
		generate_relative_pose6D_impl< use_all_pairs, false >( itBeginHand, itEndHand, std::back_inserter( dualB ) );
	}
	
	{ // bring all rotatoions into range < 180° degree
		std::transform( dualA.begin(), dualA.end(), dualA.begin(), Ubitrack::Math::hemisphere_alignment< pose_type, true >() );
		// std::transform( dualB.begin(), dualB.end(), dualB.begin(), Ubitrack::Math::substitute_rotation() );
	}
	

	const std::size_t n_cluster( std::min( n, n_select ) );
	
	std::vector< pose_type > means;
	means.reserve( n_cluster );
	Math::Stochastic::copy_probability( dualA.begin(), dualA.end(), n_cluster, std::back_inserter( means ), Math::RotationDistance< pose_type >() );
	
	
	//std::cout << "dist " << std::distance( means.begin(), means.end() );
	std::vector< std::size_t > indices;
	indices.reserve( n );
	Math::Stochastic::k_means( dualA.begin(), dualA.end(), means.begin(), means.end(), std::back_inserter( indices ), Math::RotationDistance< pose_type >() );
	
	// std::vector< pose_type > selected_a;
	// selected_a.reserve( n_cluster );
	
	// std::vector< pose_type > selected_b;
	// selected_b.reserve( n_cluster );
	
	
	typename std::vector< pose_type >::iterator itCodebook = means.begin();
	for( std::size_t i = 0; i<n_cluster; ++i, ++itCodebook )
	{
		T max_dist( 10000 );  // <- important: should be bigger than PI
		pose_type nearestPoseA;
		pose_type nearestPoseB;
		
		std::vector< size_t >::iterator itIndex = indices.begin();
		std::vector< size_t >::iterator itIndexEnd = indices.end();
		
		typename std::vector< pose_type >::iterator itPoseA = dualA.begin();
		typename std::vector< pose_type >::iterator itPoseB = dualB.begin();
		
		for( ; itIndex != itIndexEnd; ++itIndex, ++itPoseA, ++itPoseB )
		{
			if( *itIndexEnd == i )
			{
				const T d = Math::RotationDistance< pose_type >()( *itCodebook, *itPoseA );
				if( d < max_dist )
				{
					max_dist = d;
					nearestPoseA = *itPoseA;
					nearestPoseB = *itPoseB;
				}
			}
		}
		*itEye++ = Math::pose_cast< output_type >()( nearestPoseA );
		*itHand++ = Math::pose_cast< output_type >()( nearestPoseB );
		// selected_a.push_back( nearestPoseA );
		// selected_b.push_back( nearestPoseB );
	}
};

} // anonymous-namespace

UBITRACK_EXPORT void select_6DPoses( const std::vector< Math::Pose >& eyes
	, const std::vector< Math::Pose >& hands
	, const std::size_t select
	, std::vector< Math::Pose >& eyesOut
	, std::vector< Math::Pose >& handsOut  )
{
	eyesOut.reserve( select );
	handsOut.reserve( select );
	
	pose6D_selection_impl< true, double >( eyes.begin(), eyes.end(), hands.begin(), hands.end(), select, std::back_inserter( eyesOut ), std::back_inserter( handsOut ) );
	
}

/// @internal implementation function 
UBITRACK_EXPORT void generate_relative_6DPoses( const std::vector< Math::Pose >& poses, std::vector< Math::Vector< double, 8 > >& relativePoses, bool direction_flag )
{
	const std::size_t n_in = std::distance( poses.begin(), poses.end() );
	assert( n_in > 2 );
	
	// do not forget about the already inserted poses.. maybe there are some
	const std::size_t n = relativePoses.size() + (n_in*(n_in-1)) / 2 ;
	
	relativePoses.reserve( n );
	
	// choose betweenn relaitve poses for "eye" or relative poses for the "hand"
	if( direction_flag )
		generate_relative_pose6D_impl< true, true >( poses.begin(), poses.end(), std::back_inserter( relativePoses ) );
	else
		generate_relative_pose6D_impl< true, false >( poses.begin(), poses.end(), std::back_inserter( relativePoses ) );
}



}} // namespace Ubitrack::Calibration

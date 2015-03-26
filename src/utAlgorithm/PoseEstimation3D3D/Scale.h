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
 * Implementation of scale estimation in Absolute Orientation
 *
 * @author Mahmoud Baha
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */

#ifndef __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_SCALE_INCLUDED__
#define __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_SCALE_INCLUDED__


#include <cmath>	// std::sqrt
#include <cassert>	// assert-macro
#include <numeric>	// std::accumulate
#include <iterator>	// std::iterator_traits

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation3D3D {


/** @internal
 * The scale function was implemented from Mahmoud Bahaa using horn's paper 
 * 'Closed-form solution of absolute orientation using unit quaternion'
 * slide number 4 - 18/5/2010
 */
template < typename InputIterator >
typename std::iterator_traits< InputIterator >::value_type::value_type estimateScale_3D3D(
	const InputIterator iBeginLeft, const InputIterator iEndLeft
	, const InputIterator iBeginRight, const InputIterator iEndRight
	, const typename std::iterator_traits< InputIterator >::value_type& leftCentroid
	, const typename std::iterator_traits< InputIterator >::value_type& rightCentroid )
{

	/// @todo <-check that vector_type is a 3-vector
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type; 
	typedef typename vector_type::value_type value_type;
	
	// only for assertion
	{
		const std::size_t n = std::distance( iBeginLeft, iEndLeft );
		assert( n > 2u );
		const std::size_t n2 = std::distance( iBeginRight, iEndRight );
		assert( n == n2 );
	}
	
	//Compute the summation
	value_type sDenominator = 0; 
	for ( InputIterator it1 ( iBeginLeft ); it1 != iEndLeft; it1++ )
	{
		const vector_type rl_ref = (*it1) - leftCentroid;
		const value_type normSquared = (rl_ref[0]*rl_ref[0])+(rl_ref[1]*rl_ref[1])+(rl_ref[2]*rl_ref[2]);
		sDenominator += normSquared ; 	
	}	
	
	value_type sNumerator = 0;
	for ( InputIterator it2 ( iBeginRight ); it2 != iEndRight; it2++ )
	{
		const vector_type rr_ref = (*it2) - rightCentroid;
		const value_type normSquared = (rr_ref[0]*rr_ref[0])+(rr_ref[1]*rr_ref[1])+(rr_ref[2]*rr_ref[2]);						
		sNumerator += normSquared ;
	}	
	// return std::sqrt( sNumerator / sDenominator );
	return std::sqrt( sDenominator / sNumerator );
}

template < typename InputIterator >
typename std::iterator_traits< InputIterator >::value_type::value_type estimateScale_3D3D(
	const InputIterator iBeginLeft, const InputIterator iEndLeft
	, const InputIterator iBeginRight, const InputIterator iEndRight )
{

	/// @todo <-check that vector_type is a 3-vector
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type; 
	typedef typename vector_type::value_type value_type;
	
	const std::size_t n = std::distance( iBeginLeft, iEndLeft );
	
	//Compute the centroids of both coordinate systems
	const vector_type leftCentroid = std::accumulate( iBeginLeft, iEndLeft, vector_type::zeros() ) / n;
	const vector_type rightCentroid = std::accumulate( iBeginRight, iEndRight, vector_type::zeros() ) / n;
	
	// call other functions that expects the centroids, see code above
	return estimateScale_3D3D( iBeginLeft, iEndLeft, iBeginRight, iEndRight, leftCentroid, rightCentroid );
}

} } } // namespace Ubitrack::Algorithm::PoseEstimation3D3D

#endif // __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_SCALE_INCLUDED__

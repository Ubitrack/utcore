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
 * Implementation of  Absolute Orientation (3D-3D-Pose-Estimation)
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */

#include "AbsoluteOrientation.h"
#ifdef HAVE_LAPACK

#include <utMath/Blas2.h>

#include <numeric> // std::accumulate
#include <iterator> // std::iterator_traits

// namespace shortcuts
#include <boost/numeric/bindings/lapack/syev.hpp>
namespace lapack = boost::numeric::bindings::lapack;

namespace Ubitrack { namespace Calibration {

namespace { // anonymous

/** @internal */
template < class T, typename ForwardIterator >
bool calculateRotation_impl ( const ForwardIterator leftBegin, const ForwardIterator leftEnd,
		Math::Quaternion& quat,
		const ForwardIterator rightBegin, const ForwardIterator rightEnd,
		const Math::Vector< T, 3 >& leftCentroid, const Math::Vector< T, 3 >& rightCentroid )
{
	typedef typename std::iterator_traits< ForwardIterator >::value_type vector_type;
	typedef typename vector_type::value_type value_type;

	// calculate the matrix M of sums of outer-products
	Math::Matrix< T, 3, 3 > M ( Math::Matrix< T, 3, 3 >::zeros() );

	ForwardIterator itBegin1 = leftBegin;
	ForwardIterator itBegin2 = rightBegin;
	for ( ; itBegin1 != leftEnd; itBegin1++, itBegin2++ )
	{
		const vector_type vec1 = (*itBegin1) - leftCentroid;
		const vector_type vec2 = (*itBegin2) - rightCentroid;
		M += Math::OuterProduct() ( vec2, vec1 );

	}

	// calculate the matrix N as linear combinations of elements of M
	// upper right suffices, since N is symmetric
	Math::Matrix< T, 4, 4 > N;
	N( 0, 0 ) =  M ( 0, 0 )+M ( 1, 1 )+M ( 2, 2 );
	N( 1, 1 ) =  M ( 0, 0 )-M ( 1, 1 )-M ( 2, 2 );
	N( 2, 2 ) = -M ( 0, 0 )+M ( 1, 1 )-M ( 2, 2 );
	N( 3, 3 ) = -M ( 0, 0 )-M ( 1, 1 )+M ( 2, 2 );
 
	N( 0, 1 ) = M ( 1, 2 )-M ( 2, 1 );
	N( 0, 2 ) = M ( 2, 0 )-M ( 0, 2 );
	N( 0, 3 ) = M ( 0, 1 )-M ( 1, 0 );
	N( 1, 2 ) = M ( 0, 1 )+M ( 1, 0 );
	N( 1, 3 ) = M ( 2, 0 )+M ( 0, 2 );
	N( 2, 3 ) = M ( 1, 2 )+M ( 2, 1 );

	// calculate eigenvalues and eigenvectors of N
	Math::Vector< T, 4 > W;
	int i = lapack::syev ( 'V', 'U', N, W, lapack::minimal_workspace() );
	if( i != 0 )
	{
		// < 0:  if INFO = -i, the i-th argument had an illegal value
		// std::cout << "Lapack problem, eigenvalue decomposition problem " << i << ".\n";
		return false;
	}
	// largest eigenvalue is always the last one when returned by lapack
	if ( W[ 3 ] <= 0.0 )
	{
		return false;
		// UBITRACK_THROW ( "Largest Eigenvalue of Matrix N is not positive" );
	}
	
	quat = Math::Quaternion ( N ( 1, 3 ), N ( 2, 3 ), N ( 3, 3 ), N ( 0, 3 ) );
	return true; // <- return everything is ok
}

/** @internal
 * The scale function was implemented from Mahmoud Bahaa using horn's paper 
 * 'Closed-form solution of absolute orientation using unit quaternion'
 * slide number 4 - 18/5/2010
 */
template < typename ForwardIterator >
typename std::iterator_traits< ForwardIterator >::value_type::value_type calculateScale_impl(
	const ForwardIterator leftBegin, const ForwardIterator leftEnd
	, const ForwardIterator rightBegin, const ForwardIterator rightEnd )
{

	/// @todo <-check that vector_type is a 3-vector
	typedef typename std::iterator_traits< ForwardIterator >::value_type vector_type; 
	typedef typename vector_type::value_type value_type;
	
	const std::size_t n = std::distance( leftBegin, leftEnd );
	assert( n > 3 );
	assert( n == std::distance( rightBegin, rightEnd ) );
	
	//Compute the centroids of both coordinate systems
	const vector_type leftCentroid = std::accumulate( leftBegin, leftEnd, vector_type::zeros() ) / n;
	const vector_type rightCentroid = std::accumulate( rightBegin, rightEnd, vector_type::zeros() ) / n;
	
	
	//Compute the summation
	value_type sDenominator = 0; 
	for ( ForwardIterator it1 ( leftBegin ); it1 != leftEnd; it1++ )
	{
		const vector_type rl_ref = (*it1) - leftCentroid;
		const value_type normSquared = (rl_ref[0]*rl_ref[0])+(rl_ref[1]*rl_ref[1])+(rl_ref[2]*rl_ref[2]);
		sDenominator += normSquared ; 	
	}	
	
	value_type sNumerator = 0;
	for ( ForwardIterator it2 ( rightBegin ); it2 != rightEnd; it2++ )
	{
		const vector_type rr_ref = (*it2) - rightCentroid;
		const value_type normSquared = (rr_ref[0]*rr_ref[0])+(rr_ref[1]*rr_ref[1])+(rr_ref[2]*rr_ref[2]);						
		sNumerator += normSquared ;
	}	
	return std::sqrt( sNumerator / sDenominator );
}



/**
 * @internal
 * Calculate the absolute orientation problem.
 *
 * @param leftBegin iterator that points at the beginning of the 3D vectors in the left coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @param leftEnd iterator that points at one after the last of the 3D vectors in the left coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @param rightBegin iterator that points at the beginning of the 3D vectors in the right coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @param rightEnd iterator that points at one after the last of the 3D vectors in the left coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @return pose that describes the transformation of the left coordinate frame into the right coordinate frame.
 * @throws Util::Exception if lapack is not available, different number of samples for left and right side are given or the matrix N only has non-positive eigenvalues.
 */
template< typename ForwardIterator >
bool calculateAbsoluteOrientation_Impl ( const ForwardIterator itBegin1, const ForwardIterator itEnd1, Math::Pose& pose, const ForwardIterator itBegin2, const ForwardIterator itEnd2 )
{
	/// @todo <-check that vector_type is a 3-vector
	typedef typename std::iterator_traits< ForwardIterator >::value_type vector_type; 
	typedef typename vector_type::value_type value_type;

	const std::size_t n = std::distance( itBegin1, itEnd1 );
	assert( n > 3 ); // <- 3D3D needs at least three pairwise distinct values
	assert( n == std::distance( itBegin2, itEnd2 ) ); // <- should contain equal amount of values

	const vector_type leftCentroid = std::accumulate( itBegin1, itEnd1, vector_type::zeros() ) / n;
	const vector_type rightCentroid = std::accumulate( itBegin2, itEnd2, vector_type::zeros() ) / n;

	Math::Quaternion quat;
	if( !calculateRotation_impl ( itBegin1, itEnd1, quat, itBegin2, itEnd2, leftCentroid, rightCentroid ) )
		return false;
	
	const vector_type translation = leftCentroid - quat*rightCentroid;
	pose = Math::Pose ( quat, translation );
	
	return true;
}
} // anonymous-namespace

Math::Pose calculateAbsoluteOrientation ( const std::vector< Math::Vector3d >& left,
										  const std::vector< Math::Vector3d >& right)
{
	Math::Pose pose;
	// attention flipped the input parameters for new order
	calculateAbsoluteOrientation_Impl ( right.begin(), right.end(), pose, left.begin(), left.end() );
	return Math::Pose( pose );
}

double estimateScale_3D3D( const std::vector< Math::Vector3d >& m_left
	, const std::vector< Math::Vector3d >& m_right )
{
	return calculateScale_impl( m_left.begin(), m_left.end(), m_right.begin(), m_right.end() );
}

float estimateScale_3D3D( const std::vector< Math::Vector3f >& m_left
	, const std::vector< Math::Vector3f >& m_right )
{
	return calculateScale_impl( m_left.begin(), m_left.end(), m_right.begin(), m_right.end() );
}


bool estimatePose6D_3D3D( const std::vector< Math::Vector3d >& points3dA
	, Math::Pose& pose
	, const std::vector< Math::Vector3d >& points3dB )
{
	return calculateAbsoluteOrientation_Impl( points3dA.begin(), points3dA.end(), pose, points3dB.begin(), points3dB.end() );
}
	
bool estimatePose6D_3D3D( const std::vector< Math::Vector3f >& points3dA
	, Math::Pose& pose, const std::vector< Math::Vector3f >& points3dB )
{
	return calculateAbsoluteOrientation_Impl( points3dA.begin(), points3dA.end(), pose, points3dB.begin(), points3dB.end() );
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA
	, Math::Matrix3x3d& mat , const std::vector< Math::Vector3d >& points3dB )
{
	const std::size_t n = points3dA.size();
	const Math::Vector3d leftCentroid = std::accumulate( points3dA.begin(), points3dA.end(), Math::Vector3d::zeros() ) / n;
	const Math::Vector3d rightCentroid = std::accumulate( points3dB.begin(), points3dB.end(), Math::Vector3d::zeros() ) / n;

	Math::Quaternion quat;
	if( !calculateRotation_impl ( points3dA.begin(), points3dA.end(), quat, points3dB.begin(), points3dB.end(), leftCentroid, rightCentroid ) )
		return false;
	mat = Math::Matrix3x3d( quat );
	return true;
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3f >& points3dA
	, Math::Matrix3x3f& mat, const std::vector< Math::Vector3f >& points3dB )
{
	const std::size_t n = points3dA.size();
	const Math::Vector3f leftCentroid = std::accumulate( points3dA.begin(), points3dA.end(), Math::Vector3f::zeros() ) / n;
	const Math::Vector3f rightCentroid = std::accumulate( points3dB.begin(), points3dB.end(), Math::Vector3f::zeros() ) / n;

	Math::Quaternion quat;
	if( !calculateRotation_impl ( points3dA.begin(), points3dA.end(), quat, points3dB.begin(), points3dB.end(), leftCentroid, rightCentroid ) )
		return false;
	mat = Math::Matrix3x3f( quat );
	return true;
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3d >& points3dA
	, Math::Quaternion& quat, const std::vector< Math::Vector3d >& points3dB )
{
	const std::size_t n = points3dA.size();
	const Math::Vector3d leftCentroid = std::accumulate( points3dA.begin(), points3dA.end(), Math::Vector3d::zeros() ) / n;
	const Math::Vector3d rightCentroid = std::accumulate( points3dB.begin(), points3dB.end(), Math::Vector3d::zeros() ) / n;

	return calculateRotation_impl ( points3dA.begin(), points3dA.end(), quat, points3dB.begin(), points3dB.end(), leftCentroid, rightCentroid );
}

bool estimateRotation_3D3D( const std::vector< Math::Vector3f >& points3dA
	, Math::Quaternion& quat, const std::vector< Math::Vector3f >& points3dB )
{
	const std::size_t n = points3dA.size();
	const Math::Vector3f leftCentroid = std::accumulate( points3dA.begin(), points3dA.end(), Math::Vector3f::zeros() ) / n;
	const Math::Vector3f rightCentroid = std::accumulate( points3dB.begin(), points3dB.end(), Math::Vector3f::zeros() ) / n;

	return calculateRotation_impl ( points3dA.begin(), points3dA.end(), quat, points3dB.begin(), points3dB.end(), leftCentroid, rightCentroid );
}

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

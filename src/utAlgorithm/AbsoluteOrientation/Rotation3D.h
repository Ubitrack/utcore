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
 * Implementation of Absolute Orientation (3D-3D-Pose-Estimation)
 *
 * @author Manuel Huber <huberma@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */

#ifndef __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_ROTATION_3D_INCLUDED__
#define __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_ROTATION_3D_INCLUDED__

#include <utMath/Util/RotationCast.h>
#include <utMath/Blas2.h> // outer_product

#include <numeric> // std::accumulate
#include <iterator> // std::iterator_traits

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/syev.hpp>
#endif

namespace Ubitrack { namespace Algorithm { namespace AbsoluteOrientation {


/** @internal function estimating the rotation using Horn's approach*/
template < typename InputIterator, typename ResultType >
bool estimateRotation_3D3D ( const InputIterator iBeginLeft, const InputIterator iEndLeft
		, ResultType& rotation 
		, const InputIterator iBeginRight, const InputIterator iEndRight
		, const typename std::iterator_traits< InputIterator >::value_type& leftCentroid
		, const typename std::iterator_traits< InputIterator >::value_type& rightCentroid )
{
#ifndef HAVE_LAPACK
	return false;
#else
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type;
	typedef typename vector_type::value_type value_type;

	const std::size_t n = std::distance( iBeginLeft, iEndLeft );
	assert( n > 2u );
	const std::size_t n2 = std::distance( iBeginRight, iEndRight );
	assert( n == n2 );
	
	// calculate the matrix M of sums of outer-products
	Math::Matrix< value_type, 3, 3 > M ( Math::Matrix< value_type, 3, 3 >::zeros() );

	InputIterator itBegin1 = iBeginLeft;
	InputIterator itBegin2 = iBeginRight;
	for ( ; itBegin1 != iEndLeft; itBegin1++, itBegin2++ )
	{
		const vector_type vec1 = (*itBegin1) - leftCentroid;
		const vector_type vec2 = (*itBegin2) - rightCentroid;
		M += Math::outer_product( vec2, vec1 );
	}

	// calculate the matrix N as linear combinations of elements of M
	// upper right suffices, since N is symmetric
	Math::Matrix< value_type, 4, 4 > N;
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
	Math::Vector< value_type, 4 > W;
	int i = boost::numeric::bindings::lapack::syev ( 'V', 'U', N, W, boost::numeric::bindings::lapack::minimal_workspace() );
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
	
	rotation = Math::Util::RotationCast< ResultType >( )( Math::Quaternion ( N ( 1, 3 ), N ( 2, 3 ), N ( 3, 3 ), N ( 0, 3 ) ) );
	return true; // <- return everything is ok
	
#endif // HAVE_LAPACK
}

/// @internal function that also calculates the centroids
template< typename InputIterator, typename ResultType >
bool estimateRotation_3D3D( const InputIterator iBeginA, const InputIterator iEndA 
	, ResultType& result
	, const InputIterator iBeginB, const InputIterator iEndB )
{
	// some initial typedefs for easier coding
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type;
	
	const std::size_t n = std::distance( iBeginA, iEndA );
	
	const vector_type leftCentroid = std::accumulate( iBeginA, iEndA, vector_type::zeros() ) / n;
	const vector_type rightCentroid = std::accumulate( iBeginB, iEndB, vector_type::zeros() ) / n;

	return estimateRotation_3D3D ( iBeginA, iEndA, result, iBeginB, iEndB, leftCentroid, rightCentroid );
}

} } } // namespace Ubitrack::Algorithm::AbsoluteOrientation

#endif // __UBITRACK_ALGORITHM_ABSOLUTE_ORIENTATION_ROTATION_3D_INCLUDED__
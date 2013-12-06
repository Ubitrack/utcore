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

#include <utMath/Matrix.h>
#include <utUtil/Exception.h>

// namespace shortcuts
namespace ublas = boost::numeric::ublas;
#include <boost/numeric/bindings/lapack/syev.hpp>
namespace lapack = boost::numeric::bindings::lapack;

namespace Ubitrack { namespace Calibration {


/** \internal */
template < class T, typename ForwardIterator >
Math::Vector< T, 3 > calculateCentroid ( const ForwardIterator& iBegin, const ForwardIterator& iEnd, const Math::Vector< T, 3 >& )
{
	Math::Vector< T, 3 > centroid ( 0, 0, 0 );
	unsigned int count = 0;

	for ( ForwardIterator it ( iBegin ); it != iEnd; it++, count++ )
	{
		centroid += ( Math::Vector< T, 3 > )(*it);
	}
	centroid /= (T)count;
	return centroid;
}

/** \internal */
template < class T, typename ForwardIterator >
Math::Quaternion calculateRotation ( const ForwardIterator& leftBegin, const ForwardIterator& leftEnd,
									 const ForwardIterator& rightBegin, const ForwardIterator& rightEnd,
									 const Math::Vector< T, 3 >& leftCentroid, const Math::Vector< T, 3 >& rightCentroid )
{

	// calculate the matrix M of sums of products
	Math::Matrix< 3, 3, T > M ( Math::Matrix< 3, 3, T >::zeros() );

	ForwardIterator leftIterator = leftBegin;
	ForwardIterator rightIterator = rightBegin;

	for ( ; leftIterator != leftEnd; leftIterator++, rightIterator++ )
	{
		if ( rightIterator == rightEnd )
		{
			UBITRACK_THROW ( "Left side contains more samples than right side" );
		}
		M += ublas::outer_prod (*leftIterator - leftCentroid, *rightIterator - rightCentroid);
	}
	if ( rightIterator != rightEnd )
	{
		UBITRACK_THROW ( "Right Side contains more samples than left side" );
	}

	// calculate the matrix N as linear combinations of elements of M
	// upper right suffices, since N is symmetric
	Math::Matrix<4, 4, T> N;
	N (0,0) =  M ( 0, 0 )+M ( 1, 1 )+M ( 2, 2 );
	N (1,1) =  M ( 0, 0 )-M ( 1, 1 )-M ( 2, 2 );
	N (2,2) = -M ( 0, 0 )+M ( 1, 1 )-M ( 2, 2 );
	N (3,3) = -M ( 0, 0 )-M ( 1, 1 )+M ( 2, 2 );

	N (0,1) = M ( 1, 2 )-M ( 2, 1 );
	N (0,2) = M ( 2, 0 )-M ( 0, 2 );
	N (0,3) = M ( 0, 1 )-M ( 1, 0 );
	N (1,2) = M ( 0, 1 )+M ( 1, 0 );
	N (1,3) = M ( 2, 0 )+M ( 0, 2 );
	N (2,3) = M ( 1, 2 )+M ( 2, 1 );

	// calculate eigenvalues and eigenvectors of N

	Math::Vector< T, 4 > W;
	lapack::syev ( 'V', 'U', N, W, lapack::minimal_workspace() );

	// largest eigenvalue is always the last one when returned by lapack
	if ( W[3] <= 0.0 )
	{
		UBITRACK_THROW ( "Largest Eigenvalue of Matrix N is not positive" );
	}

	Math::Quaternion q ( N ( 1, 3 ), N ( 2, 3 ), N ( 3, 3 ), N ( 0, 3 ) );

	return q;
}

/** \internal */
template < class T, typename ForwardIterator >
Math::Scalar< double > calculateScale (const ForwardIterator& leftBegin, const ForwardIterator& leftEnd,
									 const ForwardIterator& rightBegin, const ForwardIterator& rightEnd
									 , const Math::Vector< T, 1 >&)
{
	//,const Math::Vector<1, T>&
	Math::Vector< T, 3 > dummy;

	//Compute the centroids of both coordinate systems
	Math::Vector< T, 3 > leftCentroid = calculateCentroid ( leftBegin, leftEnd, dummy );
	Math::Vector< T, 3 > rightCentroid = calculateCentroid ( rightBegin, rightEnd, dummy );
	
	//Reference measurements to the centroids
	Math::Vector< double, 3 > rl_ref (0.0, 0.0, 0.0);
	Math::Vector< double, 3 > rr_ref (0.0, 0.0, 0.0);
		
	double normSquared = 0.0;
	double sNumerator = 0.0;
	double sDenominator = 0.0; 
	double doubleScale = 0.0; 
	Math::Scalar< double > scale;
	
	//Compute the summation
	for ( ForwardIterator it1 ( leftBegin ); it1 != leftEnd; it1++ )
	{
		rl_ref=(Math::Vector< T, 3 >)(*it1) - leftCentroid;
		normSquared = (rl_ref[0]*rl_ref[0])+(rl_ref[1]*rl_ref[1])+(rl_ref[2]*rl_ref[2]);
		sDenominator += normSquared ; 	
	}	
	
	for ( ForwardIterator it2 ( rightBegin ); it2 != rightEnd; it2++ )
	{
		rr_ref=(Math::Vector< T, 3 >)(*it2) - rightCentroid;
		normSquared = (rr_ref[0]*rr_ref[0])+(rr_ref[1]*rr_ref[1])+(rr_ref[2]*rr_ref[2]);						
		sNumerator += normSquared ;
	}	

	doubleScale = sqrt(sNumerator/sDenominator);
	scale = Math::Scalar< double >(doubleScale);

	
	return scale;
}



/**
 * \internal
 * Calculate the absolute orientation problem.
 *
 * @param leftBegin iterator that points at the beginning of the 3D vectors in the left coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @param leftEnd iterator that points at one after the last of the 3D vectors in the left coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @param rightBegin iterator that points at the beginning of the 3D vectors in the right coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @param rightEnd iterator that points at one after the last of the 3D vectors in the left coordinate frame. Must be model of ForwardIterator and Math::Vector< double, 3 >*
 * @return pose that describes the transformation of the left coordinate frame into the right coordinate frame.
 * @throws Util::Exception if lapack is not available, different number of samples for left and right side are given or the matrix N only has non-positive eigenvalues.
 */
template < typename ForwardIterator >
Math::Pose calculateAbsoluteOrientationImpl ( const ForwardIterator leftBegin, const ForwardIterator leftEnd, const ForwardIterator rightBegin, const ForwardIterator rightEnd )
{
	Math::Vector< double, 3 > dummy;

	Math::Vector< double, 3 > leftCentroid = calculateCentroid ( leftBegin, leftEnd, dummy );
	Math::Vector< double, 3 > rightCentroid = calculateCentroid ( rightBegin, rightEnd, dummy );

	Math::Quaternion q = calculateRotation ( leftBegin, leftEnd, rightBegin, rightEnd, leftCentroid, rightCentroid );
	Math::Vector< double, 3 > translation = rightCentroid - q*leftCentroid;
	
	Math::Pose result ( q, translation );

	return result;
}
template < typename ForwardIterator >
Math::Scalar< double > calculateAbsoluteOrientationScaleImpl ( const ForwardIterator leftBegin, const ForwardIterator leftEnd, const ForwardIterator rightBegin, const ForwardIterator rightEnd )
{
	Math::Vector< double, 1 > dummyScale;

	Math::Scalar< double > scale = calculateScale( leftBegin, leftEnd, rightBegin, rightEnd, dummyScale);		  
	return scale;
}

/*
 *'calculateScale' is implemented by Mahmoud Bahaa using horn's paper 
 *'Closed-form solution of absolute orientation using unit quaternion'
 *slide number 4 - 18/5/2010
*/

Math::Pose calculateAbsoluteOrientation ( const Math::Vector< double, 3 >* leftBegin, const Math::Vector< double, 3 >* leftEnd,
										  const Math::Vector< double, 3 >* rightBegin, const Math::Vector< double, 3 >* rightEnd )
{
	return calculateAbsoluteOrientationImpl ( leftBegin, leftEnd, rightBegin, rightEnd );
}

Math::Pose calculateAbsoluteOrientation ( const std::vector< Math::Vector< double, 3 > >::iterator& leftBegin,
										  const std::vector< Math::Vector< double, 3 > >::iterator& leftEnd,
										  const std::vector< Math::Vector< double, 3 > >::iterator& rightBegin,
										  const std::vector< Math::Vector< double, 3 > >::iterator& rightEnd)
{
	return calculateAbsoluteOrientationImpl ( leftBegin, leftEnd, rightBegin, rightEnd );
}

Math::Pose calculateAbsoluteOrientation ( const std::vector< Math::Vector< double, 3 > >& left,
										  const std::vector< Math::Vector< double, 3 > >& right)
{
	return calculateAbsoluteOrientationImpl ( left.begin(), left.end(), right.begin(), right.end() );
}

Math::Scalar< double > calculateAbsoluteOrientationScale(const std::vector< Math::Vector< double, 3 > >& m_left,
															const std::vector< Math::Vector< double, 3 > >& m_right)
{
	return calculateAbsoluteOrientationScaleImpl( m_left.begin(), m_left.end(), m_right.begin(), m_right.end() );
}

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

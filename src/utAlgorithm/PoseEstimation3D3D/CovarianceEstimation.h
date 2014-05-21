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
 * Covariance estimation for absolute orientation problem
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_COVARIANCE_ESTIMATION_H_INCLUDED__
#define __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_COVARIANCE_ESTIMATION_H_INCLUDED__

// Ubitrack
#include <utCore.h>
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Blas1.h>
#include "ErrorEstimation.h"

#include <utMath/Geometry/PointTransformation.h>
#include "../Function/QuaternionRotationError.h"
#include <utMath/Stochastic/BackwardPropagation.h>

// std
#include <utility> // std::pair
#include <numeric> // std::accumulate
#include <algorithm> // std::transform

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation3D3D {


namespace Function {

template< class InputIterator >
class MultiplePointTransformationError
{
public:
	typedef typename std::iterator_traits< InputIterator >::value_type::value_type VType;
	
protected:
	const InputIterator m_iBegin;
	const InputIterator m_iEnd;
	const std::size_t m_size;
	
public:
	
	/**
	 * Constructor
	 */
	MultiplePointTransformationError( const InputIterator iBegin, const InputIterator iEnd )
		: m_iBegin( iBegin )
		, m_iEnd( iEnd )
		, m_size( std::distance( iBegin, iEnd ) )
	{}

	/**
	 * return the size of the result vector
	 */
	std::size_t size() const
	{ return 3 * m_size; }

	/**
	 * @param param containing the parameters (pose as 7-vector>
	 * @param j matrix to store the jacobian (evaluated for param) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		using namespace Ubitrack::Math;
		namespace ublas = boost::numeric::ublas;

		// convert quaternion to matrix (for speedup)
		Quaternion rotQ( Quaternion::fromVector( ublas::subrange( input, 3, 7 ) ) );
		// Matrix< VType, 3, 3 > rot( rotQ );
	
		// Matrix< VType, 3, 3 > transJ;
		Matrix< VType, 3, 3 > rotJ;
		// Vector< VType, 3 > rotated;
		// Vector< VType, 3 > translated;
		
		
		std::size_t i( 0 );
		for ( InputIterator it ( m_iBegin ); it < m_iEnd; ++it, ++i )
		{

			// // // rotate & translate points
			// noalias( rotated ) = ublas::prod( rot, *it ) + ublas::subrange( input, 0, 3 );
			// noalias( translated ) = rotated;//ublas::prod( ublas::subrange( cam, 0, 3, 0, 3 ), rotated ) + ublas::column( cam, 3 );

			// // create jacobian for this measurement
			Algorithm::Function::QuaternionRotationError< VType > qrf( *it );
			qrf.jacobian( ublas::subrange( input, 3, 7 ), rotJ );
			/// @todo _ check if jacobian is correct, needs some more review here
			Math::Matrix< VType, 3, 3 > trafo = Math::Matrix< VType, 3, 3 >::identity();
			noalias( ublas::subrange( J, i * 3, ( i + 1 ) * 3, 0, 3 ) ) 
				= trafo;
			noalias( ublas::subrange( J, i * 3, ( i + 1 ) * 3, 3, 6 ) ) 
				= ( rotJ );
		}
	}
};

} // namespace Ubitrack::Algorithm::PoseEstimation3D3D::Function



/// @internal function to calculate the covariance to an estimated pose from 3D point correspondences
template< typename InputIterator >
bool estimatePose6DCovariance( const InputIterator iBeginA, const InputIterator iEndA
	, const Math::Pose& pose
	, const InputIterator iBeginB, const InputIterator iEndB
	// , Math::Matrix< Math::Pose::value_type, 6, 6 > &covariance )
	, Math::Matrix< typename std::iterator_traits< InputIterator >::value_type::value_type, 6, 6 > &covariance )
{
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type;
	typedef typename vector_type::value_type T;

	// const std::size_t n = std::distance( iBeginA, iEndA );
	
	// const Math::Matrix< T, 3, 4 > mat( pose );

	// estimate residual error
	const T err = estimatePose6DResidual< T >( iBeginA, iEndA, pose, iBeginB, iEndB );
	
	// copy rot & trans to parameter vector
	Math::Vector< T, 7 > params;
	pose.toVector( params );

	// calculate error
	Function::MultiplePointTransformationError< InputIterator > trafoFunc( iBeginB, iEndB );

	Math::Stochastic::backwardPropagationIdentity( covariance, err, trafoFunc, params );

	return false;
}
	
}}} // namespace Ubitrack::Algorithm::PoseEstimation3D3D


#endif //__UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_COVARIANCE_ESTIMATION_H_INCLUDED__

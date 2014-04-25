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
 * Implements functions for tooltip/hotspot calibration.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */ 

// Ubitrack
#include "TipCalibration.h"
#include <utMath/Matrix.h>
#include <utMath/Blas1.h>

// std
#include <numeric>

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>
#include <boost/numeric/bindings/lapack/gels.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
# endif // HAVE_LAPACK

namespace Ubitrack { namespace Algorithm { namespace ToolTip {

#ifndef HAVE_LAPACK // write an empty body
template< typename T >
bool estimatePosition3D_6DImpl( Math::Vector< T, 3 >& 
	, const std::vector< Math::Pose >& 
	, Math::Vector< T, 3 >& )
	{ return false; }
#else // HAVE_LAPACK
/// @internal general implementation of tooltip calibration algorithm
template< typename T >
bool estimatePosition3D_6DImpl( Math::Vector< T, 3 >& pw
	, const std::vector< Math::Pose >& poses
	, Math::Vector< T, 3 >& pm )
{
	const std::size_t nPoses = ( poses.size() );
	typename Math::Matrix< T >::base_type a( 3 * nPoses, 6 );
	typename Math::Vector< T >::base_type v( 3 * nPoses );
	
	for ( std::size_t i( 0 ); i < nPoses; i++ )
	{
		// set a
		ublas::matrix_range< typename Math::Matrix< T >::base_type > r( 
			a, ublas::range( i * 3, (i+1) * 3 ), ublas::range( 0, 3 ) );
		poses[ i ].rotation().toMatrix( r );
		ublas::subrange( a, i * 3, (i+1) * 3, 3, 6 ) = - Math::Matrix< T, 3, 3 >::identity();

		// set v
		ublas::subrange( v, i * 3, (i+1) * 3 ) = -poses[ i ].translation();
	}

	// solve
	if( 0 != lapack::gels( 'N', a, v ) )
		return false;
		
	// save result
	pm = ublas::subrange( v, 0, 3 );
	pw = ublas::subrange( v, 3, 6 );
	
	return true;
}
#endif // HAVAE_LAPACK

/// @internal old function signature that performs the calibration
void tipCalibration( const std::vector< Math::Pose >& poses, 
	Math::Vector< double, 3 >& pm, Math::Vector< double, 3 >& pw )
{
	estimatePosition3D_6DImpl< double >( pw, poses, pm );
}

/// @internal specialization of tooltip calibration for type \c float
bool estimatePosition3D_6D( Math::Vector< float, 3 >& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector< float, 3 >& pm )
{
	return estimatePosition3D_6DImpl< float >( pw, poses, pm );
}

/// @internal specialization of tooltip calibration for type \c double
bool estimatePosition3D_6D( Math::Vector< double, 3 >& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector< double, 3 >& pm )
{
	return estimatePosition3D_6DImpl< double >( pw, poses, pm );
}

}}} // namespace Ubitrack::Algorithm::ToolTip

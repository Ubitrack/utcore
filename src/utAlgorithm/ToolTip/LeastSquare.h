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
 * Implements a least-square solution for tooltip/hotspot calibration.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */ 

 
#ifndef __UBITRACK_ALGROITHM_TOOLTIP_LEASTSQUARES_H_INCLUDED__
#define __UBITRACK_ALGROITHM_TOOLTIP_LEASTSQUARES_H_INCLUDED__
 
// Ubitrack
#include <utMath/Pose.h>
#include <utMath/Matrix.h>
#include <utMath/Blas1.h>

#ifdef HAVE_LAPACK

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>
#include <boost/numeric/bindings/lapack/gels.hpp>

# endif // HAVE_LAPACK

namespace Ubitrack { namespace Algorithm { namespace ToolTip {

/**
 * @ingroup tracking_algorithms
 * Computes the tooltip/hotspot calibration in a least-square fashion.
 *
 * The routine solves the following equation system using a least-square solution
 * , given a list of body @f$ i @f$ poses (R_i, t_i):
 *  @f$ (R_i -I) (p_m p_w) = -t_i @f$
 *
 * A description of the algorithm was published by Tuceryan et al. in their article
 * "Calibration requirements and procedures for a monitor-based augmented reality system"
 * in 1995 ( @cite tuceryan1995calibration ).
 *
 * @verbatim
@article{tuceryan1995calibration,
  title={Calibration requirements and procedures for a monitor-based augmented reality system},
  author={Tuceryan, Mihran and Greer, Douglas S. and Whitaker, Ross T. and Breen, David E. and Crampton, Chris and Rose, Eric and Ahlers, Klaus H},
  journal={Visualization and Computer Graphics, IEEE Transactions on},
  volume={1},
  number={3},
  pages={255--273},
  year={1995},
  publisher={IEEE}
} @endverbatim
 *
 *
 * @param pm returns the constant point in body coordinates
 * @param iBgein iterator pointing to a first pose from a container of pose elements
 * @param iEnd iterator pointing to the end of a container of pose elements
 * @param pw returns the constant point in world coordinates
 */	
template< typename T, typename InputIterator >
bool estimatePosition3D_6D( Math::Vector< T, 3 >& pw
	, const InputIterator iBegin
	, const InputIterator iEnd
	, Math::Vector< T, 3 >& pm )
{

#ifndef HAVE_LAPACK
	return false;
#else

	// shortcuts to namespaces
	namespace ublas = boost::numeric::ublas;
	namespace lapack = boost::numeric::bindings::lapack;
	
	const std::size_t nPoses = std::distance( iBegin, iEnd );
	assert( nPoses > 2 );
	
	typename Math::Matrix< T >::base_type a( 3 * nPoses, 6 );
	typename Math::Vector< T >::base_type v( 3 * nPoses );
	
	std::size_t i = 0;
	for ( InputIterator it( iBegin ); it < iEnd; ++it, ++i )
	{
		// set a
		ublas::matrix_range< typename Math::Matrix< T >::base_type > r( 
			a, ublas::range( i * 3, (i+1) * 3 ), ublas::range( 0, 3 ) );
		it->rotation().toMatrix( r );
		ublas::subrange( a, i * 3, (i+1) * 3, 3, 6 ) = - Math::Matrix< T, 3, 3 >::identity();

		// set v
		ublas::subrange( v, i * 3, (i+1) * 3 ) = -(it->translation());
	}

	// solve
	if( 0 != lapack::gels( 'N', a, v ) )
		return false;
		
	// save result
	pm = ublas::subrange( v, 0, 3 );
	pw = ublas::subrange( v, 3, 6 );
	
	return true;
#endif // HAVE_LAPACK

}

}}} // namespace Ubitrack::Algorithm::ToolTip

#endif //__UBITRACK_ALGROITHM_TOOLTIP_LEASTSQUARES_H_INCLUDED__

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
 * Implements functions for tip/hotspot calibration.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#include "TipCalibration.h"

#include <utMath/Matrix.h>

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>
#include <boost/numeric/bindings/lapack/gels.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

namespace Ubitrack { namespace Calibration {

void tipCalibration( const std::vector< Math::Pose >& poses, 
	Math::Vector< 3 >& pm, Math::Vector< 3 >& pw )
{
	const std::size_t nPoses = ( poses.size() );
	Math::Matrix< 0, 0, double >::base_type a( 3 * nPoses, 6 );
	Math::Vector< 0, double >::base_type v( 3 * nPoses );
	for ( std::size_t i( 0 ); i < nPoses; i++ )
	{
		// set a
		ublas::matrix_range< Math::Matrix< 0, 0, double >::base_type > r( 
			a, ublas::range( i * 3, (i+1) * 3 ), ublas::range( 0, 3 ) );
		poses[ i ].rotation().toMatrix( r );
		ublas::subrange( a, i * 3, (i+1) * 3, 3, 6 ) = -ublas::identity_matrix< double >( 3 );

		// set v
		ublas::subrange( v, i * 3, (i+1) * 3 ) = -poses[ i ].translation();
	}

	// solve
	lapack::gels( 'N', a, v );

	// save result
	pm = ublas::subrange( v, 0, 3 );
	pw = ublas::subrange( v, 3, 6 );
}

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

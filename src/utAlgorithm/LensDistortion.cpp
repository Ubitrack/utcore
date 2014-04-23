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
 * Implements functions for lens distortion.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#include "LensDistortion.h"
#include "Function/RadialDistortion.h"

#ifdef HAVE_LAPACK
#include <utMath/Optimization/LevenbergMarquardt.h>
#endif

namespace Ubitrack { namespace Algorithm {

template< typename T >
Math::Vector< T, 2 > projectWithDistortionImpl( const Math::Vector< T, 3 >& p, const Math::Vector< T, 4 >& dist,
	const Math::Matrix< T, 3, 3 >& K )
{
	// dehomogenize point
	Math::Vector< T, 2 > dehomogenized( p( 0 ) / p( 2 ), p( 1 ) / p( 2 ) );

	// distort
	Math::Vector< T, 2 > dp;
	Function::radialDistortion( dp, dehomogenized, dist );

	// back to image coordinates
	return Math::Vector< T, 2 >( 
		( dp( 0 ) * K( 0, 0 ) + dp( 1 ) * K( 0, 1 ) + K( 0, 2 ) ) / K( 2, 2 ), 
		(                       dp( 1 ) * K( 1, 1 ) + K( 1, 2 ) ) / K( 2, 2 ) );
}


Math::Vector< float, 2 > projectWithDistortion( const Math::Vector< float, 3 >& p, const Math::Vector< float, 4 >& dist,
	const Math::Matrix< float, 3, 3 >& K )
{
	return projectWithDistortionImpl( p, dist, K );
}


Math::Vector< double, 2 > projectWithDistortion( const Math::Vector< double, 3 >& p, const Math::Vector< double, 4 >& dist,
	const Math::Matrix< double, 3, 3 >& K )
{
	return projectWithDistortionImpl( p, dist, K );
}


template< typename T >
Math::Vector< T, 2 > lensDistortImpl( const Math::Vector< T, 2 >& p, const Math::Vector< T, 4 >& dist,
	const Math::Matrix< T, 3, 3 >& K )
{
	// unproject point to normalized camera coordinates (assuming K( 2, 2 ) == +/-1)
	T x2 = ( p( 1 ) - K( 1, 2 ) * K( 2, 2 ) ) / K( 1, 1 );
	T x1 = ( p( 0 ) - K( 0, 1 ) * x2 - K( 0, 2 ) * K( 2, 2 ) ) / K( 0, 0 );
	Math::Vector< T, 2 > camPoint( x1, x2 );

	// distort
	Math::Vector< T, 2 > dp;
	Function::radialDistortion( dp, camPoint, dist );

	// back to image coordinates
	return Math::Vector< T, 2 >( dp( 0 ) * K( 0, 0 ) + dp( 1 ) * K( 0, 1 ) + K( 0, 2 ) * K( 2, 2 ), dp( 1 ) * K( 1, 1 ) + K( 1, 2 ) * K( 2, 2 ) );
}


Math::Vector< float, 2 > lensDistort( const Math::Vector< float, 2 >& p, const Math::Vector< float, 4 >& dist,
	const Math::Matrix< float, 3, 3 >& K )
{ 
	return lensDistortImpl( p, dist, K ); 
}


Math::Vector< double, 2 > lensDistort( const Math::Vector< double, 2 >& p, const Math::Vector< double, 4 >& dist,
	const Math::Matrix< double, 3, 3 >& K )
{ 
	return lensDistortImpl( p, dist, K ); 
}

#ifdef HAVE_LAPACK

template< typename T >
Math::Vector< T, 2 > lensUnDistortImpl( const Math::Vector< T, 2 >& p, const Math::Vector< T, 4 >& dist,
	const Math::Matrix< T, 3, 3 >& K )
{
	// unproject point to normalized camera coordinates (assuming K( 2, 2 ) == +/-1)
	T x2 = ( p( 1 ) - K( 1, 2 ) * K( 2, 2 ) ) / K( 1, 1 );
	T x1 = ( p( 0 ) - K( 0, 1 ) * x2 - K( 0, 2 ) * K( 2, 2 ) ) / K( 0, 0 );
	Math::Vector< T, 2 > camPoint( x1, x2 );

	// non-linear minimization
	Math::Vector< T, 2 > dp( camPoint );
	Function::RadialDistortionWrtP< T > distFunc( dist );
	Math::Optimization::levenbergMarquardt( distFunc, dp, camPoint, Math::Optimization::OptTerminate( 5, 1e-5 ), Math::Optimization::OptNoNormalize() );

	// back to image coordinates
	return Math::Vector< T, 2 >( dp( 0 ) * K( 0, 0 ) + dp( 1 ) * K( 0, 1 ) + K( 0, 2 ) * K( 2, 2 ), dp( 1 ) * K( 1, 1 ) + K( 1, 2 ) * K( 2, 2 ) );
}


Math::Vector< float, 2 > lensUnDistort( const Math::Vector< float, 2 >& p, const Math::Vector< float, 4 >& dist,
	const Math::Matrix< float, 3, 3 >& K )
{ 
	return lensUnDistortImpl( p, dist, K ); 
}


Math::Vector< double, 2 > lensUnDistort( const Math::Vector< double, 2 >& p, const Math::Vector< double, 4 >& dist,
	const Math::Matrix< double, 3, 3 >& K )
{ 
	return lensUnDistortImpl( p, dist, K ); 
}

#endif

} } // namespace Ubitrack::Algorithm


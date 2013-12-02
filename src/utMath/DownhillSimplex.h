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
 * @ingroup math
 * @file
 * Downhill Simplex Optimizer
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 
 
#ifndef __UBITRACK_MATH_DOWNHILLSIMPLEX_INCLUDED__
#define __UBITRACK_MATH_DOWNHILLSIMPLEX_INCLUDED__

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "Optimization.h"

namespace Ubitrack { namespace Math {

/**
 * Downhill simplex minimizer á la Nelder and Mead (1965).
 * Adapted to boost::ublas from numerical recipes.
 *
 * For a better discussion of the parameters, see the levenberg-marquardt optimizer
 *
 * @param problem problem class, must have an evaluate(...) method.
 * @param params initial parameters, replaced by result at the end
 * @param measurement goal of the function
 * @param terminationCriteria when to terminate
 * @param normalize normalization function
 */
template< class P, class VT1, class VT2, class TC, class NT > 
typename VT1::value_type downhillSimplex( const P& problem, VT1& params, const VT2& measurement, 
	const TC& terminationCriteria, const NT& normalize = OptNoNormalize() )
{
	namespace ublas = boost::numeric::ublas;
	typedef typename VT1::value_type T;

	// initialize the starting points of the simplex
	unsigned ndim = params.size();
	Math::Matrix< 0, 0, VType > p( ndim + 1, ndim );
	ublas::row( p, 0 ) = params;
	for ( unsigned i = 1; i <= ndim; i++ )
	{
		ublas::matrix_row< Math::Matrix< 0, 0, T > > subrow( p, i );
		subrow = params;
		p( i, i-1 ) *= T( 1.48529 );
		normalize.evaluate( subrow, subrow );
	}

	// evaluate the function for these points
	Math::Vector< 0, T > y( ndim + 1 );
	Math::Vector< 0, T > eval( problem.size() );
	for ( unsigned i = 0; i <= ndim; i++ )
	{
		problem.evaluate( eval, ublas::row( p, i ) );
		y( i ) = ublas::norm_2( eval - measurement );
	}

	// number of iterations
	unsigned nfunk = 0;

	// the sum of all simplex points
	Math::Vector< 0, T > psum( ublas::row( p, 0 ) );
	for ( unsigned i = 1; i <= ndim; i++ )
		noalias( psum ) += ublas::row( p, i );

	while ( true )
	{
		// find lowest, highest and next-highest points
		unsigned ilo = 0;
		unsigned ihi, inhi;
		if ( y( 0 ) > y( 1 ) )
			ihi = 0, inhi = 1;
		else
			ihi = 1, inhi = 0;
		for ( unsigned i = 0; i <= ndim; i++ )
		{
			if ( y( i ) <= y( ilo ) )
				ilo = i;
			if ( y( i ) > y( ihi ) )
			{
				inhi = ihi;
				ihi = i;
			}
			else if ( y( i ) > y( inhi ) && i != ihi )
				inhi = i;
		}

		// check for termination
		if ( terminationCriteria( nfunk, y( ilo ), y( ihi ) ) )
		{
			noalias( params ) = ublas::row( p, ilo );
			return y( ilo );
		}

		nfunk += 2;

		// try extrapolation by -1.0
		T ytry = downhillSimplexTry( p, y, psum, problem, ihi, T( -1 ), eval, measurement, normalize );
		if ( ytry <= y( ilo ) )
			// better => try extrapolation by 2
			ytry = downhillSimplexTry( p, y, psum, problem, ihi, T( 2 ), eval, measurement, normalize );
		else if ( ytry >= y( inhi ) )
		{
			// worse => try intermediate lower point
			T ysave = y( ihi );
			ytry = downhillSimplexTry( p, y, psum, problem, ihi, T( 0.5 ), eval, measurement, normalize );
			if ( ytry >= ysave )
			{
				// contract around lowest point
				for ( unsigned i = 0; i <= ndim; i++ )
					if ( i != ilo )
					{
						ublas::matrix_row< Math::Matrix< 0, 0, T > > subrow( p, i );
						noalias( subrow ) += ublas::row( p, ilo );
						subrow *= T( 0.5 );
						normalize.evaluate( subrow, subrow );
						problem.evaluate( eval, subrow );
						y( i ) = ublas::norm_2( eval - measurement );
					}

				nfunk += ndim;

				// recompute psum
				noalias( psum ) = ublas::row( p, 0 );
				for ( unsigned i = 1; i <= ndim; i++ )
					noalias( psum ) += ublas::row( p, i );
			}
		}
		else
			nfunk--;
	}
}


/**
 * function internally used by the simplex optimizer
 * \internal
 */
template<  class T, class P, class VT2, class NT > 
T downhillSimplexTry( Math::Matrix< 0, 0, T >& p, Math::Vector< 0, T >& y, 
	Math::Vector< 0, T >& psum, P& problem, unsigned ihi, T fac, Math::Vector< 0, T >& eval, 
	const VT2& measurement, const NT& normalize )
{
	namespace ublas = boost::numeric::ublas;
	unsigned ndim = psum.size();
	T fac1 = ( T( 1 ) - fac ) / ndim;
	T fac2 = fac1 - fac;

	// compute the trial point
	Math::Vector< 0, T > ptry( psum * fac1 - ublas::row( p, ihi ) * fac2 );
	normalize.evaluate( ptry, ptry );

	// evaluate
	problem.evaluate( eval, ptry );
	T ytry = ublas::norm_2( eval - measurement );

	// replace the highest if the trial is better
	if ( ytry < y( ihi ) )
	{
		y( ihi ) = ytry;
		noalias( psum ) += ptry - ublas::row( p, ihi );
		noalias( ublas::row( p, ihi ) ) = ptry;
	}

	return ytry;
}


} } // namespace Ubitrack::Math

#endif // __UBITRACK_MATH_DOWNHILLSIMPLEX_INCLUDED__

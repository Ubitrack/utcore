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
 * Gauss-Newton optimizer
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_GAUSSNEWTON_INCLUDED__
#define __UBITRACK_MATH_GAUSSNEWTON_INCLUDED__

#include <utMath/Optimization/Optimization.h>

#include <utMath/Matrix.h>
#include <utMath/Vector.h>
#include <utUtil/Exception.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/bindings/blas/blas.hpp>
#include <boost/numeric/bindings/lapack/posv.hpp>
#include <boost/numeric/bindings/lapack/gels.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>


#ifdef HAVE_LAPACK

namespace Ubitrack { namespace Math {

/**
 * @ingroup math
 * Run a number of Gauss-Newton optimizer iterations.
 *
 * @par The problem class
 * The problem class P must be modeled after the UnaryFunctionPrototype and implement the function
 * \c evaluateWithJacobian which computes the predicted measurement and the jacobian wrt. the parameters to optimize.
 *
 * @param problem the problem to optimize -- provides measurement estimates and jacobians
 * @param params initial parameters on entry, optimized parameters on exit
 * @param measurement the measurement vector
 * @param normalize a UnaryFunction called after each iteration to normalize the result. Only needs to implement \c evaluate()
 * @param nIterations number of iterations
 * @return the residual of the optimization process
 */
template< class P, class VT1, class VT2, class NT >
void gaussNewton( P& problem, VT1& params, const VT2& measurement, unsigned nIterations, const NT& normalize = OptNoNormalize() )
{
	namespace lapack = boost::numeric::bindings::lapack;
	namespace blas = boost::numeric::bindings::blas;
	namespace ublas = boost::numeric::ublas;
	typedef typename VT1::value_type T;
	typedef typename Math::Matrix< T, 0, 0 > MatType;

	// create some matrices and vectors
	MatType matJacobian( measurement.size(), params.size() );
	MatType matJacobiSquare( params.size(), params.size() );
	Math::Vector< T > measurementDiff( measurement.size() );
	Math::Vector< T > paramDiff( params.size() );
	Math::Vector< T > estimatedMeasurement( measurement.size() );

	OPT_LOG_DEBUG( "Gauss-Newton entry params: " << params );

	for ( unsigned i = 0; i < nIterations; i++ )
	{
		// compute initial error
		problem.evaluateWithJacobian( estimatedMeasurement, params, matJacobian );

		ublas::noalias( measurementDiff ) = measurement - estimatedMeasurement;
		T fRes = ublas::inner_prod( measurementDiff, measurementDiff );
		OPT_LOG_TRACE( "measurementDiff: " << measurementDiff );
		OPT_LOG_DEBUG( "Gauss-Newton residual " << i << ": " << fRes );

		// start optimization loop
		blas::gemm( 'T', 'N', T( 1 ), matJacobian, matJacobian, T( 0 ), matJacobiSquare );
		blas::gemm( 'T', 'N', T( 1 ), matJacobian, measurementDiff, T( 0 ), paramDiff );

		OPT_LOG_TRACE( "Jacobian[0:16]: " << ublas::subrange( matJacobian, 0, std::min( std::size_t( 16 ), measurement.size() ), 0, params.size() ) );
		
		// do least squares
		int info = lapack::gels( 'N', matJacobiSquare, paramDiff ); // result in matParamDiff
		if ( info != 0 )
			UBITRACK_THROW( "lapack::gels returned an error" );

		noalias( params ) += paramDiff; //ublas::subrange( measurementDiff, 0, params.size() );

		OPT_LOG_TRACE( "ParamDiff (after): " << paramDiff );
		OPT_LOG_TRACE( "new params: " << params );

		// normalize
		normalize.evaluate( params, params );
	}
}

} } // namespace Ubitrack::Math

#endif // HAVE_LAPACK


#endif

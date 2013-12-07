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
 * Levenberg-Marquardt optimizer
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 
 
#ifndef __UBITRACK_MATH_LEVENBERGMARQUARDT_INCLUDED__
#define __UBITRACK_MATH_LEVENBERGMARQUARDT_INCLUDED__

#include <utMath/Optimization.h>

#ifdef HAVE_LAPACK

#include <utMath/Matrix.h>
#include <utMath/Vector.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <utUtil/Exception.h>
#include <boost/shared_ptr.hpp>

#include <boost/numeric/bindings/blas/blas.hpp>
#include <boost/numeric/bindings/lapack/posv.hpp>
#include <boost/numeric/bindings/lapack/gels.hpp>
#include <boost/numeric/bindings/lapack/gelss.hpp> // TODO: gels with blocking svd + optimal workspace query
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>



namespace Ubitrack { namespace Math {

/** possible solvers to use in levenberg-marquardt optimization */
enum LmSolverType { lmUseCholesky, lmUseQR, lmUseSVD };

/**
 * @ingroup math
 * Optimize a given problem using the levenberg marquardt optimizer.
 *
 * @par The problem class
 * The problem class P must be modeled after the UnaryFunctionPrototype and implement the function
 * \c evaluateWithJacobian which computes the predicted measurement and the jacobian wrt. the parameters to optimize.
 *
 * @param problem the problem to optimize -- provides measurement estimates and jacobians
 * @param params initial parameters on entry, optimized parameters on exit
 * @param measurement the measurement vector
 * @param terminationCriteria functor that returns true if the optimization should terminate. Is called with
 *   bool operator()( unsigned iteration, double currentError, double previousError )
 * @param normalize a UnaryFunction called after each iteration to normalize the result. Only needs to implement \c evaluate()
 * @param solver least-squares solver to use
 * @return the residual of the optimization process
 */
template< class P, class X, class Y, class TC, class NT, class WFT > 
typename X::value_type weightedLevenbergMarquardt( P& problem, X& params, const Y& measurement, 
	const TC& terminationCriteria, const NT& normalize = OptNoNormalize(), 
	 const WFT& weightFunction = OptNoWeightFunction(), LmSolverType solver = lmUseCholesky )
{
	namespace lapack = boost::numeric::bindings::lapack;
	namespace blas = boost::numeric::bindings::blas;
	namespace ublas = boost::numeric::ublas;
	typedef typename X::value_type T;
	typedef typename Math::Matrix< T, 0, 0 >::base_type MatType;
	typedef typename Math::Vector< T >::base_type VecType;
	
	// create some matrices and vectors
	boost::shared_ptr< MatType > pJacobian( new MatType( measurement.size(), params.size() ) );
	boost::shared_ptr< MatType > pJacobian2( new MatType( measurement.size(), params.size() ) );
	MatType matJacobiSquare( params.size(), params.size() );
	boost::shared_ptr< VecType > pMeasurementDiff( new VecType( measurement.size() ) );
	boost::shared_ptr< VecType > pMeasurementDiff2( new VecType( measurement.size() ) );
	VecType paramDiff( params.size() );
	VecType estimatedMeasurement( measurement.size() );
	VecType newParams( params.size() );

	// compute initial error
	problem.evaluateWithJacobian( estimatedMeasurement, params, *pJacobian );
	ublas::noalias( *pMeasurementDiff ) = measurement - estimatedMeasurement;
	OPT_LOG_TRACE( "Measurement Diff = " << *pMeasurementDiff );

	// multiply jacobian and difference with sqare root of weight matrix
	if ( !weightFunction.noWeights() )
	{
		VecType weightVector( measurement.size() );
		weightFunction.computeWeights( *pMeasurementDiff, weightVector );
		for ( unsigned i = 0; i < measurement.size(); i++ )
		{
			T w = sqrt( weightVector( i ) );
			(*pMeasurementDiff)( i ) *= w;
			ublas::row( *pJacobian, i ) *= w;
		}
		OPT_LOG_TRACE( "weights = " << weightVector );
	}

	T fErrPrev = ublas::inner_prod( *pMeasurementDiff, *pMeasurementDiff );
	OPT_LOG_DEBUG( "Levenberg-Marquardt residual 0: " << fErrPrev );

	// start optimization loop
	T fLambda = T( 1 );
	int iteration = 0;
	bool bTerminate = false;
	while ( !bTerminate )
	{
		iteration++;

		// do one optimization step
		if ( solver == lmUseCholesky )
			blas::syrk( 'L', 'T', T( 1 ), *pJacobian, T( 0 ), matJacobiSquare );
		else
			blas::gemm( 'T', 'N', T( 1 ), *pJacobian, *pJacobian, T( 0 ), matJacobiSquare );

		blas::gemm( 'T', 'N', T( 1 ), *pJacobian, *pMeasurementDiff, T( 0 ), paramDiff );
		
		// add lambda to diagonal
		for ( unsigned i = 0; i < params.size(); i++ )
			matJacobiSquare( i, i ) += fLambda;		
		
		// do least squares
		switch ( solver )
		{
		case lmUseCholesky:
			if ( lapack::posv( 'L', matJacobiSquare, paramDiff ) != 0 ) // result in paramDiff
			{
				OPT_LOG_DEBUG( "Error in cholesky decomposition, switching to SVD" );
				solver = lmUseSVD;
				continue;
			}
			break;

		case lmUseQR:
			if ( lapack::gels( 'N', matJacobiSquare, paramDiff ) != 0 ) // result in paramDiff
				UBITRACK_THROW( "lapack::gels returned an error" );
			break;

		case lmUseSVD:
			{
				Math::Vector< T > sv( params.size() );
				int rank;
				if ( lapack::gelss( matJacobiSquare, paramDiff, sv, T( -1 ), rank ) != 0 ) // result in paramDiff
					UBITRACK_THROW( "lapack::gelss returned an error" );
				OPT_LOG_DEBUG( "Effective rank: " << rank );
				OPT_LOG_TRACE( "Singular values: " << sv );
				OPT_LOG_TRACE( "Highest singular vector: " << ublas::row( matJacobiSquare, 0 ) );
				OPT_LOG_TRACE( "Lowest effective singular vector: " << ublas::row( matJacobiSquare, rank - 1 ) );
			}
			break;
		}

		OPT_LOG_TRACE( "paramDiff: " << paramDiff );
		ublas::noalias( newParams ) = params + paramDiff;

		// normalize
		normalize.evaluate( newParams, newParams );

		// compute new error
		problem.evaluateWithJacobian( estimatedMeasurement, newParams, *pJacobian2 );
		ublas::noalias( *pMeasurementDiff2 ) = measurement - estimatedMeasurement;

		// multiply jacobian and difference with sqare root of weight matrix
		if ( !weightFunction.noWeights() )
		{
			VecType weightVector( measurement.size() );
			weightFunction.computeWeights( *pMeasurementDiff2, weightVector );
			for ( unsigned i = 0; i < measurement.size(); i++ )
			{
				T w = sqrt( weightVector( i ) );
				(*pMeasurementDiff2)( i ) *= w;
				ublas::row( *pJacobian2, i ) *= w;
			}
			OPT_LOG_TRACE( "weights = " << weightVector );
		}

		T fErr = ublas::inner_prod( *pMeasurementDiff2, *pMeasurementDiff2 );

		OPT_LOG_TRACE( "measurementDiff: " << *pMeasurementDiff2 );
		OPT_LOG_DEBUG( "Levenberg-Marquardt residual " << iteration << ": " << fErr );

		// check if we should terminate
		bTerminate = terminationCriteria( iteration, fErr, fErrPrev );

		// update parameters
		if ( fErr >= fErrPrev )
			fLambda *= T( 10 );
		else
		{
			fLambda /= T( 10 );
			params = newParams;

			// swap measurementDiff
			boost::shared_ptr< VecType > pMDTemp( pMeasurementDiff );
			pMeasurementDiff = pMeasurementDiff2;
			pMeasurementDiff2 = pMDTemp;

			// swap jacobian
			boost::shared_ptr< MatType > pJTemp( pJacobian );
			pJacobian = pJacobian2;
			pJacobian2 = pJTemp;

			fErrPrev = fErr;
		}
	}

	return fErrPrev;
}

/**
 * @ingroup math
 * Optimize a given problem using the levenberg marquardt optimizer.
 *
 * @par The problem class
 * The problem class P must be modeled after the UnaryFunctionPrototype and implement the function
 * \c evaluateWithJacobian which computes the predicted measurement and the jacobian wrt. the parameters to optimize.
 *
 * @param problem the problem to optimize -- provides measurement estimates and jacobians
 * @param params initial parameters on entry, optimized parameters on exit
 * @param measurement the measurement vector
 * @param terminationCriteria functor that returns true if the optimization should terminate. Is called with
 *   bool operator()( unsigned iteration, double currentError, double previousError )
 * @param normalize a UnaryFunction called after each iteration to normalize the result. Only needs to implement \c evaluate()
 * @param solver least-squares solver to use
 * @return the residual of the optimization process
 */
template< class P, class X, class Y, class TC, class NT > 
typename X::value_type levenbergMarquardt( P& problem, X& params, const Y& measurement, 
	const TC& terminationCriteria, const NT& normalize = OptNoNormalize(), 
	LmSolverType solver = lmUseCholesky )
{ return weightedLevenbergMarquardt( problem, params, measurement, terminationCriteria, normalize, OptNoWeightFunction(), solver ); }

} } // namespace Ubitrack::Math

#endif // HAVE_LAPACK


#endif

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
 * @file
 * A vector with an associated covariance matrix.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_ERROR_COVARIANCETRANSFORM_H_INCLUDED__
#define __UBITRACK_ERROR_COVARIANCETRANSFORM_H_INCLUDED__
 
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/blas/blas3.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/ErrorVector.h>

namespace Ubitrack { namespace Math { namespace Stochastic {

/**
 * Transforms the vector by some function f and returns the result as a new vector with associated
 * covariance.
 *
 * The function to apply must be given as a function object of class F that overloads operator()
 * in the following way:
 * @verbatim
 * template< class VT1, class VT2, class MT > 
 *   void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
 * @endverbatim
 * The arguments have the following meaning:
 * - \c result: resulting vector (out)
 * - \c input: the value of the input vector (in)
 * - \c matrix object to store the jacobian in (out)
 * It can be assumed that \c VT1, \c VT2 and \c MT have the usual Boost uBlas matrix/vector semantics
 *
 * @param F function class as described above
 * @param M size of the resulting vector
 * @param f functor object of class F
 * @param resultVec vector where the result is to be stored
 * @param resultCov covariance of result
 * @param inVec input vector to be transformed
 * @param inCov covariance of input vector
 * @return a new \c ErrorVector containing the transformed vector and covariance
 */
template< class F, class VT1, class MT1, class VT2, class MT2 > 
void transformWithCovariance( const F& f, VT1& resultVec, MT1& resultCov, const VT2& inVec, const MT2& inCov )
{
	namespace ublas = boost::numeric::ublas;
	typedef typename VT1::value_type VType;

	// compute result and jacobian
	Math::Matrix< VType, 0, 0 > jacobian( resultVec.size(), inVec.size() );
	f.evaluateWithJacobian( resultVec, inVec, jacobian );
	
	// transform covariance
	Math::Matrix< VType, 0, 0 > im( ublas::prod( jacobian, inCov ) );
	noalias( resultCov ) = ublas::prod( im, ublas::trans( jacobian ) );
}


/** 
 * overload for \c ErrorVector objects.
 * Note: the size of the result (M) must be given explicitly as template parameter!
 */
template< unsigned M, unsigned N, class F, typename VType  > inline
ErrorVector< VType, M > transformWithCovariance( const F& f, const ErrorVector< VType, N > in )
{
	ErrorVector< VType, M > result;
	transformWithCovariance( f, result.value, result.covariance, in.value, in.covariance );
	return result;
}


/**
 * Apply a function f that updates a subvector of a given vector from another subvector 
 * of the same vector and update the covariance matrix.
 * 
 * Example usage: time update of a kalman filter.
 *
 * Note: The subvectors may overlap.
 *
 * @param f, see \c transform
 * @param value the vector to be transformed
 * @param covariance covariance of the vector
 * @param iOutBegin first index of the subvector to be updated
 * @param iOutEnd index of element after the subvector to be updated
 * @param iInBegin first index of subvector used as input to f
 * @param iInEnd index of element after subvector used as input to f
 */
template< class F, class VT, class MT > 
void transformRangeInternalWithCovariance( const F& f, VT& value, MT& covariance, 
	unsigned iOutBegin, unsigned iOutEnd, unsigned iInBegin, unsigned iInEnd )
{
	namespace ublas = boost::numeric::ublas;
	typedef typename VT::value_type VType;
	
	const unsigned nOutSize = iOutEnd - iOutBegin;
	const unsigned nInSize = iInEnd - iInBegin;
	Math::Vector< VType > result( nOutSize );
	
	// compute result and jacobian
	Math::Matrix< VType, 0, 0 > jacobian( nOutSize, nInSize );
	f.evaluateWithJacobian( result, ublas::subrange( value, iInBegin, iInEnd ), jacobian );
	ublas::subrange( value, iOutBegin, iOutEnd ) = result;
	
	// transform covariance in a block-matrix fashion
	Math::Matrix< VType, 0, 0 > im( ublas::prod( jacobian, ublas::subrange( covariance, iInBegin, iInEnd, 0, value.size() ) ) );
	
	// the left/upper row/column
	if ( iOutBegin > 0 )
	{
		ublas::subrange( covariance, iOutBegin, iOutEnd, 0, iOutBegin ) = 
			ublas::subrange( im, 0, nOutSize, 0, iOutBegin );
		noalias( ublas::subrange( covariance, 0, iOutBegin, iOutBegin, iOutEnd ) ) = 
			ublas::trans( ublas::subrange( im, 0, nOutSize, 0, iOutBegin ) );
	}
	
	// the diagonal part
	noalias( ublas::subrange( covariance, iOutBegin, iOutEnd, iOutBegin, iOutEnd ) ) = 
		ublas::prod( ublas::subrange( im, 0, nOutSize, iInBegin, iInEnd ), ublas::trans( jacobian ) );
		
	// the right/lower row/column
	if ( iOutEnd < value.size() )
	{
		ublas::subrange( covariance, iOutBegin, iOutEnd, iOutEnd, value.size() ) = 
			ublas::subrange( im, 0, nOutSize, iOutEnd, value.size() );
		noalias( ublas::subrange( covariance, iOutEnd, value.size(), iOutBegin, iOutEnd ) ) = 
			ublas::trans( ublas::subrange( im, 0, nOutSize, iOutEnd, value.size() ) );
	}		
}


/** Overload for \c ErrorVector */
template< unsigned N, class F, class VType > inline 
void transformRangeInternalWithCovariance( const F& f, ErrorVector< VType, N >& v,
	unsigned iOutBegin, unsigned iOutEnd, unsigned iInBegin, unsigned iInEnd )
{
	transformRangeInternalWithCovariance( f, v.value, v.covariance, iOutBegin, iOutEnd, iInBegin, iInEnd );
}


/**
 * Transforms two vectors by some function f and computes the result with covariance
 *
 */
template< class F, class VT1, class MT1, class VT2, class MT2, class VT3, class MT3 > 
void binarytransformWithCovariance( const F& f, VT1& resultVec, MT1& resultCov, 
	const VT2& inVec1, const MT2& inCov1, const VT3& inVec2, const MT3& inCov2 )
{
	namespace ublas = boost::numeric::ublas;
	typedef typename VT1::value_type VType;

	// compute result and jacobian
	Math::Matrix< VType, 0, 0 > jacobian1( resultVec.size(), inVec1.size() );
	Math::Matrix< VType, 0, 0 > jacobian2( resultVec.size(), inVec2.size() );
	f.evaluateWithJacobian( resultVec, inVec1, inVec2, jacobian1, jacobian2 );
	
	// transform covariance
	Math::Matrix< VType, 0, 0 > im1( ublas::prod( jacobian1, inCov1 ) );
	noalias( resultCov ) = ublas::prod( im1, ublas::trans( jacobian1 ) );
	Math::Matrix< VType, 0, 0 > im2( ublas::prod( jacobian2, inCov2 ) );
	noalias( resultCov ) += ublas::prod( im2, ublas::trans( jacobian2 ) );
}


/** overload for \c ErrorVector objects */
template< unsigned M, unsigned N, unsigned K, class F, class VType  > inline
ErrorVector< VType, M > binarytransformWithCovariance( const F& f, const ErrorVector< VType, N > in1, const ErrorVector< VType, K > in2 )
{
	ErrorVector< VType, M > result;
	binarytransformWithCovariance( f, result.value, result.covariance, in1.value, in1.covariance, in2.value, in2.covariance );
	return result;
}


}}} // namespace Ubitrack::Math::Stochastic

#endif

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
 * Backward propagation of covariance
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 
 
#ifndef __UBITRACK_MATH_BACKWARDPROPAGATION_INCLUDED__
#define __UBITRACK_MATH_BACKWARDPROPAGATION_INCLUDED__



#ifdef HAVE_LAPACK

#include <utMath/Matrix.h>
#include <utMath/Vector.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/bindings/blas/blas.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>

namespace Ubitrack { namespace Math { namespace Stochastic {


/**
 * @ingroup math
 * Performs backward propagation of covariance.
 * 
 * FIXME: this function doesn't seeem to work correctly. 
 * Use backwardPropagationIdentity, which gives different (and more plausible) results.
 *
 * Given a function y = f(x), known parameters x and a covariance matrix E of y, the backward propagation
 * computes the covariance C of the parameters x. size(E) >= size(C).
 *
 * @param result matrix C where the resulting covariance is stored
 * @param input matrix containing the input covariance E
 * @param function class modeled after \c UnaryFunctionPrototype, describes the measurement function f
 * @param params parameters x of the function f.
 */
template< class MT1, class MT2, class F, class VT1 > 
void backwardPropagation( MT1& result, const MT2& input, const F& function, const VT1& params )
{
	// !!!! FIXME: this function doesn't seeem to work.  !!!!!
	// Use backwardPropagationIdentity, which gives different results
	
	namespace lapack = boost::numeric::bindings::lapack;
	namespace blas = boost::numeric::bindings::blas;
	namespace ublas = boost::numeric::ublas;

	typedef typename MT1::value_type VType;
	typedef typename Math::Matrix< VType, 0, 0 > MatType;

	// evaluate jacobian
	MatType jacobian( input.size1(), result.size1() );
	function.jacobian( params, jacobian );
	
	// factorize E into Q^T * T * Q
	MatType Q( input );
	Math::Vector< VType > t( input.size1() );
	lapack::syev( 'V', 'U', Q, t, lapack::minimal_workspace() );
	
	// compute Q * J
	MatType QJ( jacobian.size1(), jacobian.size2() );
	blas::gemm( Q, jacobian, QJ );

    // compute T^(-1/2)	
	VType precision( t( 0 ) * VType( 1e-8 ) );
	for ( unsigned i = 0; i < result.size1(); i++ )
		if ( t( i ) < precision )
			t( i ) = 0;
		else
			t( i ) = 1 / sqrt( t( i ) );

	// compute J' = T^(-1/2) * Q * J
	for ( unsigned i = 0; i < result.size1(); i++ )
		ublas::row( QJ, i ) *= t( i );

	// perform SVD on QJ
	MatType dummyU( function.size(), function.size() );
	MatType dummyVt( result.size1(), result.size1() );
	Math::Vector< VType > s( result.size1() );
	lapack::gesvd( 'N', 'O', QJ, s, dummyU, dummyVt );

	// compute singular values of pseudo-inverse
	precision = s( 0 ) * VType( 1e-8 );
	for ( unsigned i = 0; i < s.size(); i++ )
		if ( s( i ) < precision )
			s( i ) = 0;
		else
			s( i ) = 1 / s( i );

	// compute covariance
	for ( unsigned i = 0; i < s.size(); i++ )
		ublas::row( QJ, i ) *= s( i );

	blas::syrk( 'U', 'T', VType( 1 ), ublas::subrange( QJ, 0, s.size(), 0, s.size() ), VType( 0 ), result );

	// fill lower half of matrix
	for ( unsigned r = 1; r < result.size1(); r++ )
		for ( unsigned c = 0; c < r; c++ )
			result( r, c ) = result( c, r );
}
		

		

/**
 * @ingroup math
 * Performs backward propagation of covariance where the input covariance is identity multiplied by some factor \c s.
 *
 * Given a function y = f(x), known parameters x and a covariance matrix E = s * I of y, the backward propagation
 * computes the covariance C of the parameters x. size(E) >= size(C).
 *
 * @param result matrix C where the resulting covariance is stored
 * @param s scaling factor of input covariance E. E = s * I
 * @param jacobian Jacobian matrix of the function f evaluated at x.
 *   Note: the matrix must be column_major and will be modified!
 */
template< class MT1, class MT3 > 
void backwardPropagationIdentity( MT1& result, typename MT1::value_type s, MT3& jacobian )
{
	namespace lapack = boost::numeric::bindings::lapack;
	namespace blas = boost::numeric::bindings::blas;
	namespace ublas = boost::numeric::ublas;

	typedef typename MT1::value_type VType;
	typedef typename Math::Matrix< VType, 0, 0 > MatType;

	// perform SVD on the jacobian
	MatType dummyU( jacobian.size1(), jacobian.size1() );
	MatType dummyVt( result.size1(), result.size1() );
	Math::Vector< VType > v( result.size1() );
	lapack::gesvd( 'N', 'O', jacobian, v, dummyU, dummyVt );

	// compute singular values of pseudo-inverse, multiplied by sqrt(s)
	VType precision( v( 0 ) * VType( 1e-8 ) );
	VType sSqrt = sqrt( s );
	for ( unsigned i = 0; i < v.size(); i++ )
		if ( v( i ) < precision )
			v( i ) = 0;
		else
			v( i ) = sSqrt / v( i );

	// compute covariance
	for ( unsigned i = 0; i < v.size(); i++ )
		ublas::row( jacobian, i ) *= v( i );

	blas::syrk( 'U', 'T', VType( 1 ), ublas::subrange( jacobian, 0, v.size(), 0, v.size() ), VType( 0 ), result );
	
	// fill lower half of matrix
	for ( unsigned r = 1; r < result.size1(); r++ )
		for ( unsigned c = 0; c < r; c++ )
			result( r, c ) = result( c, r );
}


/**
 * @ingroup math
 * Performs backward propagation of covariance where the input covariance is identity multiplied by some factor \c s.
 *
 * Given a function y = f(x), known parameters x and a covariance matrix E = s * I of y, the backward propagation
 * computes the covariance C of the parameters x. size(E) >= size(C).
 *
 * @param result matrix C where the resulting covariance is stored
 * @param s scaling factor of input covariance E. E = s * I
 * @param function class modeled after \c UnaryFunctionPrototype, describes the measurement function f
 * @param params parameters x of the function f.
 */
template< class MT1, class F, class VT1 > 
void backwardPropagationIdentity( MT1& result, typename MT1::value_type s, const F& function, const VT1& params )
{
	namespace ublas = boost::numeric::ublas;

	typedef typename MT1::value_type VType;
	typedef typename Math::Matrix< VType, 0, 0 > MatType;

	// evaluate jacobian J
	MatType jacobian( function.size(), result.size1() );
	function.jacobian( params, jacobian );
	
	backwardPropagationIdentity( result, s, jacobian );
}


/**
 * @ingroup math
 * Performs backward propagation of covariance where the input covariance E is a diagonal matrix given by
 * the vector \c e.
 *
 * Given a function y = f(x), known parameters x and a covariance matrix E = diag( e ) of y, the backward propagation
 * computes the covariance C of the parameters x. size(E) >= size(C).
 *
 * @param result matrix C where the resulting covariance is stored
 * @param input vector e containing the diagonal elements of the input covariance E
 * @param jacobian of the function f evaluated at x. Note: must be column_major and will be modified.
 */
template< class MT1, class VT2, class MT3 > 
void backwardPropagationDiagonal( MT1& result, const VT2& input, MT3& jacobian )
{
	namespace lapack = boost::numeric::bindings::lapack;
	namespace blas = boost::numeric::bindings::blas;
	namespace ublas = boost::numeric::ublas;

	typedef typename MT1::value_type VType;
	typedef typename Math::Matrix< VType, 0, 0 > MatType;

	// multiply jacobian with E^(-1/2)
	for ( unsigned i = 0; i < result.size1(); i++ )
		ublas::row( jacobian, i ) *= VType( 1 ) / sqrt( input( i ) );

	// perform SVD on jacobian
	MatType dummyU( jacobian.size1(), jacobian.size1() );
	MatType dummyVt( result.size1(), result.size1() );
	Math::Vector< VType > s( result.size1() );
	lapack::gesvd( 'N', 'O', jacobian, s, dummyU, dummyVt );

	// compute singular values of pseudo-inverse
	VType precision = s( 0 ) * VType( 1e-8 );
	for ( unsigned i = 0; i < s.size(); i++ )
		if ( s( i ) < precision )
			s( i ) = 0;
		else
			s( i ) = 1 / s( i );

	// compute covariance
	for ( unsigned i = 0; i < s.size(); i++ )
		ublas::row( jacobian, i ) *= s( i );

	blas::syrk( 'U', 'T', VType( 1 ), ublas::subrange( jacobian, 0, s.size(), 0, s.size() ), VType( 0 ), result );

	// fill lower half of matrix
	for ( unsigned r = 1; r < result.size1(); r++ )
		for ( unsigned c = 0; c < r; c++ )
			result( r, c ) = result( c, r );
}

		
/**
 * @ingroup math
 * Performs backward propagation of covariance where the input covariance E is a diagonal matrix given by
 * the vector \c e.
 *
 * Given a function y = f(x), known parameters x and a covariance matrix E = diag( e ) of y, the backward propagation
 * computes the covariance C of the parameters x. size(E) >= size(C).
 *
 * @param result matrix C where the resulting covariance is stored
 * @param input vector e containing the diagonal elements of the input covariance E
 * @param function class modeled after \c UnaryFunctionPrototype, describes the measurement function f
 * @param params parameters x of the function f.
 */
template< class MT1, class VT2, class F, class VT1 > 
void backwardPropagationDiagonal( MT1& result, const VT2& input, const F& function, const VT1& params )
{
	namespace ublas = boost::numeric::ublas;

	typedef typename MT1::value_type VType;
	typedef typename Math::Matrix< VType, 0, 0 > MatType;

	// evaluate jacobian
	MatType jacobian( input.size(), result.size1() );
	function.jacobian( params, jacobian );

	backwardPropagationDiagonal( result, input, jacobian );
}
		
}}} // namespace Ubitrack::Math::Stochastic

#endif // HAVE_LAPACK

#endif

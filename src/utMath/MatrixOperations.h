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
 * Inverse functions for boost::numeric::ublas::matrix.
 * @author Daniel Muhra <muhra@in.tum.de>
 */


#ifndef _Ubitrack_Math_MatrixOperations_INCLUDED_
#define _Ubitrack_Math_MatrixOperations_INCLUDED_

#include "Vector.h"
#include "Matrix.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>


#include <utUtil/StaticAssert.h>

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#endif


namespace Ubitrack { namespace Math {

#ifdef HAVE_LAPACK

/**
 * compute the determinant of a matrix m
 * @param T a matrix type
 * @param m a matrix of type T
 * @return the matrix determinant
 */
template< class T > 
typename T::value_type determinant( const T& mat )
{
	namespace ublas = boost::numeric::ublas;

	// make a copy of mat, as the factorization will overwrite the contents	
	Math::Matrix< 0, 0, typename T::value_type > a( mat );
	Math::Vector< 0, int > ipiv( a.size1() );
	typedef typename T::size_type size_type;
	
	boost::numeric::bindings::lapack::getrf( a, ipiv );
	
	typename T::value_type det = 1;

    for ( size_type i=0; i < ipiv.size(); ++i) 
	{
		if ( size_type( ipiv( i ) ) != i+1 )
			det *= -1;
		det *= a( i, i );
    }
    return det;
}

/**
 * invert the square matrix m
 * @param T a matrix type
 * @param m a matrix of type T
 * @return the inverted matrix
 */
template< class T > T invert_matrix( const T& m )
{
	// make a copy of m, as the factorization will overwrite the contents
	T a( m );
	Math::Vector< 0, int > ipiv( a.size1() );

	// factorize and compute inverse
	boost::numeric::bindings::lapack::getrf( a, ipiv );
	boost::numeric::bindings::lapack::getri( a, ipiv );

	return a;
}

/**
 * compute the pseudo-inverse  of the matrix m
 * @param M a matrix type
 * @param mat a matrix of type M
 * @return the inverted matrix
 */
template< typename T, typename ST > 
boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, ST > pseudoInvert_matrix( const boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, ST >& mat )
{
	namespace ublas = boost::numeric::ublas;
	typedef ublas::matrix< T, ublas::column_major, ST > M;
	typedef typename M::size_type size_type;
	// make a copy of m, as the factorization will overwrite the contents
	M a( mat );
	
	size_type n = mat.size1();
	size_type m = mat.size2();
	size_type nSingularValues = std::min( n, m );
	
	Math::Vector< 0, T > s( nSingularValues );
	Math::Matrix< 0, 0, T > U( n, n );
	Math::Matrix< 0, 0, T > Vt( m, m );

	boost::numeric::bindings::lapack::gesvd( 'S', 'S', a, s, U, Vt );

	for( size_type i=0; i<nSingularValues; i++ )
	{
		if( s( i ) != 0 )
		{
			s( i ) = static_cast< T >( 1 ) / s( i );
		}
	}

	if( n > m )
	{
		for ( size_type i = 0; i < m; i++ )
			ublas::row( Vt, i ) *= s( i );

		a = ublas::prod( ublas::subrange( U, 0, n, 0, m ) , Vt );
	}
	else
	{
		for ( size_type i = 0; i < n; i++ )
			ublas::column( U, i ) *= s( i );

		a = ublas::prod( U , ublas::subrange( Vt, 0, n, 0, m ) );
	}

	return M( ublas::trans( a ) );
}

#endif // HAVE_LAPACK

} } // namespace Ubitrack::Math

#endif // _Ubitrack_Math_Matrix_INCLUDED_


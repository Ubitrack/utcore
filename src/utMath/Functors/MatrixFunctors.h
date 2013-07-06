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
 * Functors for common matrix operations
 *
 * The Functors can easily be applied to containers like 
 * std::vector< Math::Matrix< N, N, T > > using std::transform.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 
 
#ifndef __H__MATRIX_FUNCTORS__
#define __H__MATRIX_FUNCTORS__


// Ubitrack
#include <utMath/Matrix.h>
#include <utMath/MatrixOperations.h>

namespace Ubitrack { namespace Math { namespace Functors {


/**
 * @ingroup math
 * Functor to calculate the determinant of a matrix.
 */
struct matrix_determinant
{
public:

	/**
	 * @ingroup math
	 * Calculates the determinant of a 2x2 matrix.
	 *
	 * @param matrix 2x2 matrix
	 * @return determinant of the 2x2 matrix
	 */
	template< typename T >
	T operator() ( const Math::Matrix< 2, 2, T >& matrix ) const
    {
		const T A1_1 = matrix( 0, 0 );
		const T A1_2 = matrix( 0, 1 );
		const T A2_1 = matrix( 1, 0 );
		const T A2_2 = matrix( 1, 1 );
		return ( matrix( 0, 0 )*matrix( 1, 1 ) - matrix( 0, 1 )*matrix( 1, 0 ) );
    }
	
	/**
	 * @ingroup math
	 * Calculates the determinant of a 3x3 matrix.
	 *
	 * @param matrix 3x3 matrix
	 * @return determinant of the 3x3 matrix
	 */
	template< typename T >
	T operator() ( const Math::Matrix< 3, 3, T >& matrix ) const
    {
		const T val1 = matrix( 0, 0 )*matrix( 1, 1 )*matrix( 1, 1 );
		const T val2 = matrix( 0, 0 )*matrix( 1, 2 )*matrix( 1, 2 );
		const T val3 = matrix( 1, 2 )*matrix( 1, 0 )*matrix( 1, 1 );
		const T val4 = matrix( 1, 2 )*matrix( 1, 2 )*matrix( 2, 0 );
		const T val5 = matrix( 0, 2 )*matrix( 1, 0 )*matrix( 1, 2 );
		const T val6 = matrix( 0, 2 )*matrix( 1, 1 )*matrix( 2, 0 );
		
		return val1-val2-val3+val4+val5-val6;
    }
	
#ifdef HAVE_LAPACK
	/**
	 * @ingroup math
	 * Calculates the determinant of a nxn matrix.
	 *
	 * @param matrix nxn matrix
	 * @return determinant of the nxn matrix
	 */
	template< std::size_t N, typename T >
	T operator() ( const Math::Matrix< N, N, T >& matrix ) const
    {
		return Math::determinant( matrix );
	}
#endif
};


/**
 * @ingroup math
 * Functor to calculate the inverse of a matrix.
 */
struct matrix_inverse
{

public:

	/**
	 * @ingroup math
	 * Calculates the inverse of a 2x2 matrix.
	 *
	 * @param matrix 2x2 matrix
	 * @return inverse of the 2x2 matrix
	 */
	template< typename T >
	Math::Matrix< 2, 2, T > operator() ( const Math::Matrix< 2, 2, T >& matrix ) const
    {
		const T A1_1 = matrix( 0, 0 );
		const T A1_2 = matrix( 0, 1 );
		const T A2_1 = matrix( 1, 0 );
		const T A2_2 = matrix( 1, 1 );
		const T determinant(A1_1*A2_2-A1_2*A2_1);
		
		T pInverse[ 4 ]; //(row-major)
		pInverse[ 0 ] =  A2_2/determinant;
		pInverse[ 1 ] = -A1_2/determinant;
		pInverse[ 2 ] = -A2_1/determinant;
		pInverse[ 3 ] =  A1_1/determinant;
		return Math::Matrix< 2, 2, T >( pInverse );
    }
	
	/**
	 * @ingroup math
	 * Calculates the inverse of a 3x3 matrix.
	 *
	 * @param matrix 3x3 matrix
	 * @return inverse of the 3x3 matrix
	 */
	template< typename T >
    const Math::Matrix< 3, 3, T > operator() ( const Math::Matrix< 3, 3, T >& matrix ) const
    {
		const T A1_1 = matrix( 0, 0 );
		const T A1_2 = matrix( 0, 1 );
		const T A1_3 = matrix( 0, 2 );
		
		const T A2_1 = matrix( 1, 0 );
		const T A2_2 = matrix( 1, 1 );
		const T A2_3 = matrix( 1, 2 );
		
		const T A3_1 = matrix( 2, 0 );
		const T A3_2 = matrix( 2, 1 );
		const T A3_3 = matrix( 2, 2 );
		
		const T determinant (A1_1*A2_2*A3_3-A1_1*A2_3*A3_2-A1_2*A2_1*A3_3+A1_2*A2_3*A3_1+A1_3*A2_1*A3_2-A1_3*A2_2*A3_1);
		
		T pInverse[ 9 ]; //(row-major)
		pInverse[ 0 ] = (A2_2*A3_3-A2_3*A3_2) / determinant;
		pInverse[ 1 ] = -(A1_2*A3_3-A1_3*A3_2) / determinant;
		pInverse[ 2 ] = (A1_2*A2_3-A1_3*A2_2) / determinant;
		pInverse[ 3 ] = -(A2_1*A3_3-A2_3*A3_1) / determinant;
		pInverse[ 4 ] = (A1_1*A3_3-A1_3*A3_1) / determinant;
		pInverse[ 5 ] = -(A1_1*A2_3-A1_3*A2_1) / determinant;
		pInverse[ 6 ] = (A2_1*A3_2-A2_2*A3_1) / determinant;
		pInverse[ 7 ] = -(A1_1*A3_2-A1_2*A3_1) / determinant;
		pInverse[ 8 ] = (A1_1*A2_2-A1_2*A2_1) / determinant;
		
		return Math::Matrix< 3, 3, T >( pInverse );
    }
#ifdef HAVE_LAPACK
	/**
	 * @ingroup math
	 * Calculates the inverse of a nxn matrix.
	 *
	 * @param matrix nxn matrix
	 * @return inverse of the nxn matrix
	 */
	template< std::size_t N, typename T >
	Math::Matrix< N, N, T > operator() ( const Math::Matrix< N, N, T >& matrix ) const
    {
		return Math::invert_matrix( matrix );
	}
#endif
};

} } } //Ubitrack::Math::Functors

#endif // __H__MATRIX_FUNCTORS__

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
 * Wrapper for boost::numeric::ublas::matrix.
 * @author Florian Echtler <echtler@in.tum.de>
 */


#ifndef _Ubitrack_Math_Matrix_INCLUDED_
#define _Ubitrack_Math_Matrix_INCLUDED_

// WARNING: all boost/serialization headers should be 
//          included AFTER all boost/archive headers
#include <boost/serialization/access.hpp>
#include <boost/numeric/ublas/matrix.hpp>


#include <utUtil/StaticAssert.h>

#include "Quaternion.h"
#include "Vector.h"
#include "Pose.h"

namespace Ubitrack { namespace Math {

	
// forward declaration of Matrix class
template< int M, int N, typename T = double > class Matrix;


/// stream output operator
template< int M, int N, typename T > std::ostream& operator<<( std::ostream& s, const Matrix< M, N, T >& m )
{
	for( int i = 0; i < M; i++ ) {
		s << "[ ";
		for( int j = 0; j < N; j++ )
			s << m(i,j) << " ";
		s << "]\n";
	}
	return s;
}

/**
 * @ingroup math
 * Wraps a boost::numeric::ublas::matrix for convenience.
 * Provides stream output, serialization and additional mathematical operations.
 * Note: we are storing the matrix column-major for compatibility with fortran/lapack. Sorry for any confusion this may cause.
 * @param M number of rows
 * @param N number of columns
 * @param T type (defaults to double)
 */
template< int M, int N, typename T > class Matrix
 	: public boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, boost::numeric::ublas::bounded_array< T, M*N > >
{
	friend std::ostream& operator<< <> ( std::ostream& s, const Matrix< M, N, T >& m );
	friend class ::boost::serialization::access;

	public:
		typedef boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, boost::numeric::ublas::bounded_array< T, M*N > > BaseMatrix;
	
		/** default constructor */
		Matrix()
			: BaseMatrix( M, N )
		{ }
		
		/** pro-forma constructor */
		Matrix( unsigned size1, unsigned size2 )
			: BaseMatrix( M, N )
		{ assert( size1 == M && size2 == N ); }

		/**
		 * construct from ublas matrix expression
		 * @param e a matrix_expression
		 */
		template< class ME > 
		Matrix( const boost::numeric::ublas::matrix_expression< ME >& e )
			: BaseMatrix( e )
		{ assert( e().size1() == M && e().size2() == N ); }

		/**
		 * construct from array 
		 * @param pFirst an array with M*N elements (row-major)
		 */
		Matrix( const T* pFirst )
			: BaseMatrix( M, N )
		{
			for( int i = 0; i < M; i++ )
				for( int j = 0; j < N; j++ )
					(*this)(i,j) = pFirst[i*N+j];
		}

		/**
		 * create a 4*4 homogeneous pose matrix
		 * @param rotation a quaternion
		 * @param position a translation vector
		 */
		Matrix( const Quaternion& rotation, const Vector< 3 >& position )
			: BaseMatrix( M, N )
		{
			UBITRACK_STATIC_ASSERT( ((M==3)||(M==4))&&(N==4), MATRIX_4x4_REQUIRED );

			rotation.toMatrix( *this );

			(*this)(0,3) = T( position[0] );
			(*this)(1,3) = T( position[1] );
			(*this)(2,3) = T( position[2] );

			if(M>3)
			{
				(*this)(3,0) = 0;
				(*this)(3,1) = 0;
				(*this)(3,2) = 0;
				(*this)(3,3) = 1;
			}
		}

		/**
		 * create a 4*4 homogeneous pose matrix
		 * @param pose a pose
		 */
		Matrix( const Pose& pose )
			: BaseMatrix( M, N )
		{
			UBITRACK_STATIC_ASSERT( ((M==4)||(M==3))&&(N==4), MATRIX_3x4_OR_4x4_REQUIRED );

			pose.rotation().toMatrix( *this );

			Vector< 3 > position = pose.translation();

			(*this)(0,3) = T( position[0] );
			(*this)(1,3) = T( position[1] );
			(*this)(2,3) = T( position[2] );

			if(M>3)
			{
				(*this)(3,0) = 0;
				(*this)(3,1) = 0;
				(*this)(3,2) = 0;
				(*this)(3,3) = 1;
			}
		}

		/**
		 * create a 3*3 rotation matrix
 		 * @param rotation a quaternion
		 */
		Matrix( const Quaternion& rotation )
			: BaseMatrix( M, N )
		{
			UBITRACK_STATIC_ASSERT( (M==3)&&(N==3), MATRIX_3x3_REQUIRED );

			rotation.toMatrix( *this );
		}

		/** might facilitate compiler optimizations */
		std::size_t size1() const
		{ return M; }
		
		/** might facilitate compiler optimizations */
		std::size_t size2() const
		{ return N; }
	
		
		/** return a pointer to the raw data */
		T* content() 
		{ return (*this).data().begin(); }

		/** return a pointer to the raw data */
		const T* content() const
		{ return (*this).data().begin(); }

		/** length of the raw data */
		int size() 
		{ return (*this).data().size();  }

		/**
		 * assign from ublas matrix_expression
		 * @param e a matrix_expression
		 */
		template< class AE > 
		Matrix< M, N, T >& operator=( const boost::numeric::ublas::matrix_expression< AE >& e )
		{
			assert( e().size1() == M && e().size2() == N );
			BaseMatrix::operator=( e );
			return *this;
		}

		/** 
		 * assign from base class
		 * @param m a boost::numeric::ublas::matrix
		 */
		Matrix< M, N, T >& operator=( const BaseMatrix& m )
		{
			assert( m().size1() == M && m.size2() == N );
			BaseMatrix::operator=( m );
			return *this;
		}

        /**
         * @return N*N square identity matrix
         */
        static Matrix< M, N, T> identity()
        {
            return  boost::numeric::ublas::identity_matrix< T >( M );
        }

        /**
         * @return M*N zero matrix
         */
        static Matrix< M, N, T> zeros()
        {
            return  boost::numeric::ublas::zero_matrix< T >( M, N );
        }

	protected:

		/// serialize this Matrix
		template< class Archive > 
		void serialize( Archive& ar, const unsigned int version )
		{
			for( int i = 0; i < M; i++ )
				for( int j = 0; j < N; j++ )
					ar & (*this)(i,j);
		}

};

/** 
 * converts a left-hand matrix into a right hand matrix
 * the old matrix will be overwritten
 * @param matrix the matrix to be converted
 */

template< class M >
void leftHandToRightHandMatrix(M& matrix)
{
	UBITRACK_STATIC_ASSERT( (M==3)&&(N==3), MATRIX_3x3_REQUIRED );

	matrix( 2, 0 ) *= -1;
	matrix( 2, 1 ) *= -1;
	matrix( 2, 3 ) *= -1;

	matrix( 0, 2 ) *= -1;
	matrix( 1, 2 ) *= -1;
	matrix( 3, 2 ) *= -1;
}

/** compares two matrices */
template< int M, int N, typename T >
bool operator==( const Matrix< M, N, T >& a, const Matrix< M, N, T >& b )
{
	for ( unsigned r = 0; r < M; r++ )
		for ( unsigned c = 0; c < N; c++ )
			if ( a( r, c ) != b( r, c ) )
				return false;
	return true;
}

// Some common used types
// Note that the type qualifier is dropped and double is used

typedef Math::Matrix < 2, 2 > Matrix2x2;
typedef Math::Matrix < 3, 3 > Matrix3x3;
typedef Math::Matrix < 4, 4 > Matrix4x4;
typedef Math::Matrix < 3, 4 > Matrix3x4;


} } // namespace Ubitrack::Math


// define traits for boost lapack bindings
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/detail/generate_const.hpp>

namespace boost { namespace numeric { namespace bindings { namespace traits {

template< int sM, int sN, typename T, typename M >
struct matrix_detail_traits< Ubitrack::Math::Matrix< sM, sN, T >, M > 
	: matrix_detail_traits< typename Ubitrack::Math::Matrix< sM, sN, T >::BaseMatrix, typename detail::generate_const< M, typename M::BaseMatrix >::type >
{
};

} } } } // namespace boost::numeric::bindings::traits

#endif // _Ubitrack_Math_Matrix_INCLUDED_


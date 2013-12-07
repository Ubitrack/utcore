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


#ifndef __MATRIX_H_INCLUDED__
#define __MATRIX_H_INCLUDED__

#include <utUtil/StaticAssert.h>

#include "Quaternion.h"
#include "Vector.h"
#include "Pose.h"

#include <cstddef> //  std::size_t
//next include is only for output stuff -> remove to own header?
#include <iostream> //std::ostream

// WARNING: all boost/serialization headers should be 
//          included AFTER all boost/archive headers
#include <boost/serialization/access.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace Ubitrack { namespace Math {

	



/**
 * @ingroup math
 * Wraps a boost::numeric::ublas::matrix for convenience.
 * Provides stream output, serialization and additional mathematical operations.
 * Note: we are storing the matrix column-major for compatibility with fortran/lapack. Sorry for any confusion this may cause.
 * @tparam M number of rows
 * @tparam N number of columns
 * @tparam T type (defaults to \c double)
 */
template< typename T, std::size_t M = 0, std::size_t N = M > class Matrix
 	: public boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, boost::numeric::ublas::bounded_array< T, M*N > >
{
	public:
		typedef boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, boost::numeric::ublas::bounded_array< T, M*N > > base_type;
		typedef Math::Matrix< T, M, N >	self_type;
		typedef T						value_type;
		typedef typename base_type::size_type size_type;
		
		/** default constructor */
		Matrix()
			: base_type( M, N )
		{ }
		
		/** pro-forma constructor */
		Matrix( size_type size1, size_type size2 )
			: base_type( M, N )
		{ assert( size1 == M && size2 == N ); }

		/**
		 * construct from ublas matrix expression
		 * @param e a matrix_expression
		 */
		template< class ME > 
		Matrix( const boost::numeric::ublas::matrix_expression< ME >& e )
			: base_type( e )
		{ assert( e().size1() == M && e().size2() == N ); }

		/**
		 * construct from array 
		 * @param pFirst an array with M*N elements (row-major)
		 */
		Matrix( const T* pFirst )
			: base_type( M, N )
		{
			for( size_type i = 0; i < M; i++ )
				for( size_type j = 0; j < N; j++ )
					(*this)(i,j) = pFirst[i*N+j];
		}

		/**
		 * create a 4*4 homogeneous pose matrix
		 * @param rotation a quaternion
		 * @param position a translation vector
		 */
		Matrix( const Quaternion& rotation, const Vector< double, 3 >& position )
			: base_type( M, N )
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
			: base_type( M, N )
		{
			UBITRACK_STATIC_ASSERT( ((M==4)||(M==3))&&(N==4), MATRIX_3x4_OR_4x4_REQUIRED );

			pose.rotation().toMatrix( *this );

			Vector< double, 3 > position = pose.translation();

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
			: base_type( M, N )
		{
			UBITRACK_STATIC_ASSERT( (M==3)&&(N==3), MATRIX_3x3_REQUIRED );

			rotation.toMatrix( *this );
		}

		/** might facilitate compiler optimizations */
		size_type size1() const
		{ return M; }
		
		/** might facilitate compiler optimizations */
		size_type size2() const
		{ return N; }
	
		
		/** return a pointer to the raw data */
		T* content() 
		{ return (*this).data().begin(); }

		/** return a pointer to the raw data */
		const T* content() const
		{ return (*this).data().begin(); }

		/** length of the raw data */
		size_type size() 
		{ return (*this).data().size();  }

		/**
		 * assign from ublas matrix_expression
		 * @param e a matrix_expression
		 */
		template< class AE > 
		Matrix< T, M, N >& operator=( const boost::numeric::ublas::matrix_expression< AE >& e )
		{
			assert( e().size1() == M && e().size2() == N );
			base_type::operator=( e );
			return *this;
		}

		/** 
		 * assign from base class
		 * @param m a boost::numeric::ublas::matrix
		 */
		Matrix< T, M, N >& operator=( const base_type& m )
		{
			assert( m().size1() == M && m.size2() == N );
			base_type::operator=( m );
			return *this;
		}

        /**
         * @return N*N square identity matrix
         */
        static Matrix< T, M, N > identity()
        {
            return  boost::numeric::ublas::identity_matrix< T >( M );
        }

        /**
         * @return M*N zero matrix
         */
        static Matrix< T, M, N > zeros()
        {
            return  boost::numeric::ublas::zero_matrix< T >( M, N );
        }

	protected:

		friend class ::boost::serialization::access;
		
		/// serialize this Matrix
		template< class Archive > 
		void serialize( Archive& ar, const unsigned int version )
		{
			for( size_type i = 0; i < M; i++ )
				for( size_type j = 0; j < N; j++ )
					ar & (*this)(i,j);
		}

};


/**
 * @ingroup math
 * Specialization of Matrix-class for matrices of varying sizes.
 * The size of the matrix is determined at runtime.
 * Many functions are dropped since they are not neccessary yet( e.g. serialize ).
 * @tparam T type of elements, defaults to \c double
 */
template< typename T >
class Matrix< T, 0, 0 >
 	: public boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, boost::numeric::ublas::unbounded_array< T > >
{
	public:
		
		typedef boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major, boost::numeric::ublas::unbounded_array< T > > base_type;
		typedef Math::Matrix< T, 0, 0 >			self_type;
		typedef T								value_type;
		typedef typename base_type::size_type	size_type;
		
		/** Default constructor */
		Matrix( )
			: base_type( )
		{ }
		
		/**
		 * Constructor from given dimensions
		 *
		 * @param size1 first dimension (=M) of the M-by-N matrix (rows)
		 * @param size2 second dimension (=N) of the M-by-N matrix (columns)
		 */
		Matrix( const size_type size1, const size_type size2 )
			: base_type( size1, size2 )
		{  }
		/**
		 * construct from ublas matrix expression
		 * @param e a matrix_expression
		 */
		template< class ME > 
		Matrix( const boost::numeric::ublas::matrix_expression< ME >& e )
			: base_type( e )
		{ }
		/**
		 * assign from ublas matrix_expression
		 * @param e a matrix_expression
		 */
		template< class AE > 
		Matrix< T, 0, 0 >& operator=( const boost::numeric::ublas::matrix_expression< AE >& e )
		{
			base_type::operator=( e );
			return *this;
		}
		/**
		 * Zero Matrix
		 *
		 * All values of the matrix will be set to zero.
		 * @param size1 first dimension (=M) of the M-by-N zero matrix (rows)
		 * @param size2 second dimension (=N) of the M-by-N zero matrix (columns)
         * @return M-by-N zero matrix
         */
        static Matrix< T, 0, 0 > zeros( const size_type size1, const size_type size2 )
        {
            return boost::numeric::ublas::zero_matrix< T >( size1, size2 );
        }
		
		/**
		 * Identity Matrix
		 *
		 * All values of the diagonal of the matrix will be set to one,
		 * all other values will be set to zero.
		 * @param size dimension (=N) of the N-by-N identity matrix.
         * @return N-by-N identity matrix
         */
        static Matrix< T, 0, 0 > identity( const size_type size )
        {
            return boost::numeric::ublas::identity_matrix< T >( size );
        }
		
		/**
		 * Scalar Matrix
		 *
		 * All values of the matrix will be set to the scalar value.
		 * @param size1 first dimension (=M) of the M-by-N scalar matrix (rows)
		 * @param size2 second dimension (=N) of the M-by-N scalar matrix (columns)
		 * @param value the scalar that is set on all matrix values
         * @return M-by-N scalar matrix
         */
		static Matrix< T, 0, 0 > scalar( const size_type size1, const size_type size2, const T value )
        {
            return boost::numeric::ublas::scalar_matrix< T >( size1, size2, value );
        }
};

/// stream output operator for a single Matrix
template<  typename T, std::size_t M, std::size_t N >
std::ostream& operator<<( std::ostream& s, const Matrix< T, M, N >& m )
{
	for( std::size_t i = 0; i < M; i++ ) {
		s << "[ ";
		for( std::size_t j = 0; j < N; j++ )
			s << m(i,j) << " ";
		s << "]\n";
	}
	return s;
}

/** 
 * converts a left-hand matrix into a right hand matrix
 * the old matrix will be overwritten
 * @param matrix the matrix to be converted
 */

template< class Ma >
void leftHandToRightHandMatrix(Ma& matrix)
{
	// TODO: Assert does not work this way..
	// UBITRACK_STATIC_ASSERT( (M==3)&&(N==3), MATRIX_3x3_REQUIRED );

	matrix( 2, 0 ) *= -1;
	matrix( 2, 1 ) *= -1;
	matrix( 2, 3 ) *= -1;

	matrix( 0, 2 ) *= -1;
	matrix( 1, 2 ) *= -1;
	matrix( 3, 2 ) *= -1;
}

/** compares two matrices */
template<  typename T, std::size_t M, std::size_t N >
bool operator==( const Matrix< T, M, N >& a, const Matrix< T, M, N >& b )
{
	for ( std::size_t r = 0; r < M; r++ )
		for ( std::size_t c = 0; c < N; c++ )
			if ( a( r, c ) != b( r, c ) )
				return false;
	return true;
}

/// typedef for a 2-by-2 Matrix of type \c double
typedef Math::Matrix< double, 2, 2 > Matrix2x2d;
/// typedef for a 3-by-3 Matrix of type \c double
typedef Math::Matrix< double, 3, 3 > Matrix3x3d;
/// typedef for a 4-by-4 (transformation) Matrix of type \c double
typedef Math::Matrix< double, 4, 4 > Matrix4x4d;
/// typedef for a 3-by-4 (projection) Matrix of type \c double
typedef Math::Matrix< double, 3, 4 > Matrix3x4d;

/// typedef for a 2-by-2 Matrix of type \c float
typedef Math::Matrix< float, 2, 2 > Matrix2x2f;
/// typedef for a 3-by-3 Matrix of type \c float
typedef Math::Matrix< float, 3, 3 > Matrix3x3f;
/// typedef for a 4-by-4 (transformation) Matrix of type \c float
typedef Math::Matrix< float, 4, 4 > Matrix4x4f;
/// typedef for a 3-by-4 (projection) Matrix of type \c float
typedef Math::Matrix< float, 3, 4 > Matrix3x4f;


} } // namespace Ubitrack::Math


// define traits for boost lapack bindings
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/detail/generate_const.hpp>

namespace boost { namespace numeric { namespace bindings { namespace traits {

template< std::size_t sM, std::size_t sN, typename T, typename M >
struct matrix_detail_traits< Ubitrack::Math::Matrix< T, sM, sN >, M > 
	: matrix_detail_traits< typename Ubitrack::Math::Matrix< T, sM, sN >::base_type, typename detail::generate_const< M, typename M::base_type >::type >
{
};

} } } } // namespace boost::numeric::bindings::traits

#endif // __MATRIX_H_INCLUDED__

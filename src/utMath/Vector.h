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
 * Wrapper for boost::numeric::ublas::vector.
 * @author Daniel Pustka <pustka@in.tum.de>
 */


#ifndef __VECTOR_H_INCLUDED__
#define __VECTOR_H_INCLUDED__

#include <utUtil/StaticAssert.h>

#include <assert.h>
#include <cstddef> //  std::size_t
//next three includes are only for output stuff -> remove to own header?
#include <iterator> //std::outstream_iterator
#include <iostream> //std::ostream
#include <algorithm> //std::copy


// WARNING: all boost/serialization headers should be
//          included AFTER all boost/archive headers
#include <boost/serialization/access.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace Ubitrack { namespace Math {

/// forward declaration of Vector class with default template parameter
template< std::size_t N, class T = double > class Vector;

/**
 * @ingroup math
 * Wraps a boost::numeric::ublas::vector for convenience.
 * Provides stream output and serialization as well as additional mathematical operations.
 * @param N size of the vector
 * @param T type of elements, defaults to double
 */
template< std::size_t N, class T > class Vector
	: public boost::numeric::ublas::vector< T, boost::numeric::ublas::bounded_array< T, N > >
{
	public:
		typedef boost::numeric::ublas::vector< T, boost::numeric::ublas::bounded_array< T, N > > base_type;
		typedef Math::Vector< N, T >	self_type;
		typedef T						value_type;
		typedef std::size_t				size_type;

		/** Default constructor */
		Vector()
			: base_type( N )
		{}
		
		/** pro-forma constructor */
		Vector( const std::size_t size )
			: base_type( N )
		{ assert( size == N ); }

		/** construct from fixed number of parameters (2) */
		Vector( T p0, T p1 )
			: base_type( N )
		{
			UBITRACK_STATIC_ASSERT( N == 2, VECTOR_2_REQUIRED );
			(*this)[0] = p0;
			(*this)[1] = p1;
		}

		/** construct from fixed number of parameters (3) */
		Vector( T p0, T p1, T p2 )
			: base_type( N )
		{
			UBITRACK_STATIC_ASSERT( N == 3, VECTOR_3_REQUIRED );
			(*this)[0] = p0;
			(*this)[1] = p1;
			(*this)[2] = p2;
		}

		/** construct from fixed number of parameters (4) */
		Vector( T p0, T p1, T p2, T p3 )
			: base_type( N )
		{
			UBITRACK_STATIC_ASSERT( N == 4, VECTOR_4_REQUIRED );
			(*this)[0] = p0;
			(*this)[1] = p1;
			(*this)[2] = p2;
			(*this)[3] = p3;
		}

		/**
		 * Construct from vector_expression
		 * @param e a vector_expression
		 */
		template< class AE > 
		Vector( const boost::numeric::ublas::vector_expression< AE >& e )
			: base_type( e )
		{ assert( e().size() == N ); }

		/**
		 * Construct from array.
		 * Note: This will copy the contents.
		 * @param pFirst pointer to first element
		 */
		Vector( const T* pFirst )
			: base_type( N )
		{
			for( size_type i = 0; i < N; i++ )
				(*this)[i]= pFirst[i];
		}

		/** might facilitate compiler optimizations */
		size_type size() const
		{ return N; }
		
		/** return a pointer to the raw data */
		T* content() 
		{ return (*this).data().begin(); }
		
		/** return a pointer to the raw data */
		const T* content() const
		{ return (*this).data().begin(); }
		
		/**
		 * assign from vector_expression
		 * @param e a vector_expression
		 */
		template< class AE > 
		Vector< N, T >& operator=( const boost::numeric::ublas::vector_expression< AE >& e )
		{
			assert( e().size() == N );
			base_type::operator=( e );
			return *this;
		}

		/**
		 * assign from base class
		 * @param v a ublas::vector
		 */
		Vector< N, T >& operator=( const base_type& v )
		{
			assert( v().size() == N );
			base_type::operator=( v );
			return *this;
		}

		Vector< N, T >& operator<( const base_type& v )
		{
			assert( v().size() == N );
			base_type::operator<( v );
			return *this;
		}

		bool operator<( const Vector< N, T >& v ) const
		{
			for( size_type i = 0; i < N; i++ )
			{
				if((*this)[i] < v[i])	return true;
				if((*this)[i] > v[i])	return false;
			}
			return false;
		}
		/**
		 * Returns vector of zeros
		 */
		static Vector< N, T> zeros()
		{
			return boost::numeric::ublas::zero_vector< T >( N );
		}
		
	protected:
	
		friend class ::boost::serialization::access;
		
		/** serialize Vector object */
		template< class Archive > 
		void serialize( Archive& ar, const unsigned int version )
		{
			for( size_type i = 0; i < N; i++ )
				ar & (*this)[i];
		}
};

/**
 * @ingroup math
 * Specialization of Vector-class for vectors of varying sizes.
 * The size of the vector is determined at runtime.
 * Many functions are dropped since they are not necessary yet ( e.g. serialize ).
 * @tparam T type of elements, defaults to double
 */
template< typename T >
class Vector< 0, T >
	: public boost::numeric::ublas::vector< T, boost::numeric::ublas::unbounded_array< T > >
{
	public:
		// some typedefs for templated algorithm design
		typedef boost::numeric::ublas::vector< T, boost::numeric::ublas::unbounded_array< T > > base_type;
		typedef Math::Vector< 0, T >			self_type;
		typedef T								value_type;
		typedef typename base_type::size_type	size_type;

		/** Default constructor */
		Vector( )
			: base_type( )
		{}

		/** Constructor from a given dimension. */
		Vector( const size_type size )
			: base_type( size )
		{ }

		/**
		 * Construct from vector_expression
		 * @param e a vector_expression
		 */
		template< class AE > 
		Vector( const boost::numeric::ublas::vector_expression< AE >& e )
			: base_type( e )
		{ }

		/**
		 * assign from vector_expression
		 * @param e a vector_expression
		 */
		template< class AE > 
		Vector< 0, T >& operator=( const boost::numeric::ublas::vector_expression< AE >& e )
		{
			base_type::operator=( e );
			return *this;
		}

		/**
         * Returns Vector of zeros
         */
        static Vector< 0, T > zeros( const size_type size )
        {
            return boost::numeric::ublas::zero_vector< T >( size );
        }
};

/** stream output operator for a single Math::Vector */
template< std::size_t N, typename T >
std::ostream& operator<<( std::ostream& s, const Vector< N, T >& v )
{
	s << "[ ";
	for( std::size_t i = 0; i < N; i++ )
		s << v[i] << " ";
	s << "]";
	return s;
}

/** stream output operator for any (stl-)container of Math::Vector */
template< std::size_t N, typename T, template < typename V_TYPE, typename V_ALLOC > class container >
std::ostream& operator<<( std::ostream& s, const container< Math::Vector< N, T >, std::allocator< Math::Vector< N, T > > >& vec_cont )
{
	s << "[ ";
	std::copy( vec_cont.begin(), vec_cont.end(), std::ostream_iterator< Math::Vector< N, T > > ( s, " " ) );
	s << "]";
	return s;
}

/**
 * computes a linear interpolation between two Vectors of same size
 * @param x first Vector
 * @param y second Vector
 * @param t interpolation point between 0.0 and 1.0
 * @return the interpolated vector
 */
template< std::size_t N, class T >
Vector < N, T > linearInterpolate ( const Vector< N, T >& x, const Vector< N, T >& y, double t ) 
{
	T w1 = 1.0 - t;
	T w2 = t;

	return Vector< N, T >( w1*x + w2*y );
}


/**
 * computes the cross product of two 3-vectors
 * @param a first Vector
 * @param b second Vector
 * @return the cross product of a and b
 */
template< class E1, class E2 >
Vector< 3, typename E1::value_type > cross_prod( 
	const boost::numeric::ublas::vector_expression< E1 >& a, const boost::numeric::ublas::vector_expression< E2 >& b )
{
	Vector< 3, typename E1::value_type > r;
	r( 0 ) = a()( 1 ) * b()( 2 ) - a()( 2 ) * b()( 1 );
	r( 1 ) = a()( 2 ) * b()( 0 ) - a()( 0 ) * b()( 2 );
	r( 2 ) = a()( 0 ) * b()( 1 ) - a()( 1 ) * b()( 0 );
	return r;
}


/** compares two vectors */
template< std::size_t N, typename T >
bool operator==( const Vector< N, T >& a, const Vector< N, T >& b )
{
	for ( std::size_t c = 0; c < N; c++ )
		if ( a( c ) != b( c ) )
			return false;
	return true;
}

/// typedef for 2 dimensional Vector of type \c double
typedef Math::Vector< 2, double > Vector2d;
/// typedef for 3 dimensional Vector of type \c double
typedef Math::Vector< 3, double > Vector3d;
/// typedef for 4 dimensional Vector of type \c double
typedef Math::Vector < 4, double > Vector4d;


/// typedef for 2 dimensional Vector of type \c float
typedef Math::Vector< 2, float > Vector2f;
/// typedef for 3 dimensional Vector of type \c float
typedef Math::Vector< 3, float > Vector3f;
/// typedef for 4 dimensional Vector of type \c float
typedef Math::Vector< 4, float > Vector4f;

} } // namespace Ubitrack::Math


// define traits for boost lapack bindings
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/detail/generate_const.hpp>

namespace boost { namespace numeric { namespace bindings { namespace traits {

template< std::size_t sN, typename T, typename M >
struct vector_detail_traits< Ubitrack::Math::Vector< sN, T >, M > 
	: vector_detail_traits< typename Ubitrack::Math::Vector< sN, T >::base_type, typename detail::generate_const< M, typename M::base_type >::type >
{
};

} } } } // namespace boost::numeric::bindings::traits

#endif // __VECTOR_H_INCLUDED__

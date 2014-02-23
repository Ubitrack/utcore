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

#ifndef __UBTRACK_MATH_VECTOR_H_INCLUDED__
#define __UBTRACK_MATH_VECTOR_H_INCLUDED__

#include <utUtil/StaticAssert.h>

#include <assert.h>
#include <cstddef> //  std::size_t
//next three includes are only for output stuff -> remove to own header?
#include <iosfwd> //std::ostream
// #include <iterator> //std::ostream_iterator
#include <algorithm> //std::copy


// WARNING: all boost/serialization headers should be
//          included AFTER all boost/archive headers
#include <boost/serialization/access.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace Ubitrack { namespace Math {

/**
 * @ingroup math
 * @brief A general \b Vector data structure for mathematical purposes.
 *
 * This template is designed to provide a mathemtical vector
 * data-structure for the various mathematical functions within Ubitrack
 * and store vectorized information.
 
 * Two template sepcializations for different types of programming are provided.
 * One can either specify the size (=dimension) of the \b Vector at compile time
 * as a second template arguement and therefore dafeine a stack allocated Vector
 * for fast computations. Alternatively defining a size at compile time can be
 * skipped, such that the \b Vector of the given dimension will be allocated 
 * dynamically on the heap.
 *
 * @note Several common mathematical functions for a Vector are provided
 * in \c Blas1.h, \c Blas2.h and \c VectorFunctions.h .
 * @note Wraps a boost::numeric::ublas::vector for convenience.
 * @note Provides stream output and serialization as well as some additional mathematical operations.
 * @attention see Vector< T, 0 > for the second specialization.
 *
 * @tparam T builtin type of Vector elements (e.g \c double or \c float )
 * @tparam N size of vector (defaults to zero for a heap allocated \b Vector)
 */
template< typename T, std::size_t N = 0 > class Vector
	: public boost::numeric::ublas::vector< T, boost::numeric::ublas::bounded_array< T, N > >
{
	public:
		typedef boost::numeric::ublas::vector< T, boost::numeric::ublas::bounded_array< T, N > > base_type;
		typedef Math::Vector< T, N >	self_type;
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
		Vector< T, N >& operator=( const boost::numeric::ublas::vector_expression< AE >& e )
		{
			assert( e().size() == N );
			base_type::operator=( e );
			return *this;
		}

		/**
		 * assign from base class
		 * @param v a ublas::vector
		 */
		Vector< T, N >& operator=( const base_type& v )
		{
			assert( v().size() == N );
			base_type::operator=( v );
			return *this;
		}

		Vector< T, N >& operator<( const base_type& v )
		{
			assert( v().size() == N );
			base_type::operator<( v );
			return *this;
		}

		bool operator<( const Vector< T, N >& v ) const
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
		static Vector< T, N > zeros()
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
 * @brief Specialization of \b Vector for memory allocation during runtime.
 *
 * The size of the \b Vector is determined at runtime, the heap is used
 * for memory allocation in this case.
 * Many functions are dropped since they are not necessary yet ( e.g. \c serialize ).
 *
 * @note Please see Vector for more details on the stack allocated version.
 *
 * @tparam T builtin type of vector elements (e.g \c double or \c float )
 */
template< typename T >
class Vector< T, 0 >
	: public boost::numeric::ublas::vector< T, boost::numeric::ublas::unbounded_array< T > >
{
	public:
		// some typedefs for templated algorithm design
		typedef boost::numeric::ublas::vector< T, boost::numeric::ublas::unbounded_array< T > > base_type;
		typedef Math::Vector< T, 0 >			self_type;
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
		Vector< T, 0 >& operator=( const boost::numeric::ublas::vector_expression< AE >& e )
		{
			base_type::operator=( e );
			return *this;
		}

		/**
         * Returns Vector of zeros
         */
        static Vector< T, 0 > zeros( const size_type size )
        {
            return boost::numeric::ublas::zero_vector< T >( size );
        }
};

/** @internal stream output operator for a single Math::Vector of a specific dimension */
template< typename T, std::size_t N >
std::ostream& operator<<( std::ostream& s, const Vector< T, N >& v )
{
	s << "[ ";
	for( std::size_t i = 0; i < N; ++i )
		s << v[ i ] << " ";
	s << "]";
	return s;
}


/** @internal specialization for stream output operator for a single Math::Vector of a (at compile time) unknown dimension */
template< typename T >
std::ostream& operator<<( std::ostream& s, const Vector< T >& v )
{
	const std::size_t n ( v.size() );
	s << "[ ";
	for( std::size_t i = 0; i<n; ++i )
		s << v[ i ] << " ";
	s << "]";
	return s;
}

/** @internal stream output operator for any (stl-)container of Math::Vector */
template< typename T, std::size_t N, template < typename , typename > class container >
std::ostream& operator<<( std::ostream& s, const container< Math::Vector< T, N >, std::allocator< Math::Vector< T, N > > >& vec_cont )
{
	s << "[ ";
	std::copy( vec_cont.begin(), vec_cont.end(), std::ostream_iterator< Math::Vector< T, N > > ( s, " " ) );
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
template< typename T, std::size_t N >
Vector < T, N > linearInterpolate ( const Vector< T, N >& x, const Vector< T, N >& y, double t ) 
{
	T w1 = 1.0 - t;
	T w2 = t;

	return Vector< T, N >( w1*x + w2*y );
}

/** compares two vectors */
template< typename T, std::size_t N >
bool operator==( const Vector< T, N >& a, const Vector< T, N >& b )
{
	for ( std::size_t c = 0; c < N; c++ )
		if ( a( c ) != b( c ) )
			return false;
	return true;
}

/// typedef for 2 dimensional Vector of type \c double
typedef Math::Vector< double, 2 > Vector2d;
/// typedef for 3 dimensional Vector of type \c double
typedef Math::Vector< double, 3 > Vector3d;
/// typedef for 4 dimensional Vector of type \c double
typedef Math::Vector < double, 4 > Vector4d;


/// typedef for 2 dimensional Vector of type \c float
typedef Math::Vector< float, 2 > Vector2f;
/// typedef for 3 dimensional Vector of type \c float
typedef Math::Vector< float, 3 > Vector3f;
/// typedef for 4 dimensional Vector of type \c float
typedef Math::Vector< float, 4 > Vector4f;

} } // namespace Ubitrack::Math


// define traits for boost lapack bindings
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/detail/generate_const.hpp>

namespace boost { namespace numeric { namespace bindings { namespace traits {

template< std::size_t sN, typename T, typename M >
struct vector_detail_traits< Ubitrack::Math::Vector< T, sN >, M > 
	: vector_detail_traits< typename Ubitrack::Math::Vector< T, sN >::base_type, typename detail::generate_const< M, typename M::base_type >::type >
{
};

} } } } // namespace boost::numeric::bindings::traits

#endif // __UBTRACK_MATH_VECTOR_H_INCLUDED__

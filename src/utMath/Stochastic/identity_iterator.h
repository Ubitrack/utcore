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
 * @ingroup util
 * @file
 * identity (iterator) template.
 * 
 * This container/iterator adds the possibility to easily
 * use stl-algorithms, like \c std::transform, with a
 * constant element.
 *
 * @todo File will be moved later to another place. -> do not include more
 * often than necessary.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __UBITRACK_UTIL_IDENTITY_ITERATOR_H__
#define __UBITRACK_UTIL_IDENTITY_ITERATOR_H__

#include <iterator> // std::iterator, std::random_access_iterator_tag

namespace Ubitrack { namespace Util { 

/**
 * @ingroup util
 * identity (iterator) template.
 * 
 * This container/iterator adds the possibility to easily
 * use stl-algorithms, like \c std::transform, with a
 * constant element.
 *
 * @tparam T the type of the element the container should store.
 */ 
template< typename T >
struct identity
{
	typedef identity< T >										self_type;
	
	T&															value;
	
	size_t														counter;
	
	identity( T& value_in )
		: value ( value_in )
		, counter( 0 )
	{}
	
	identity( T& value_in, const std::size_t n_in )
		: value ( value_in )
		, counter( n_in )
	{ }
	
	/** @todo change this to random access iterator or maybe reverse_iterator */
	struct iterator
		: public std::iterator< std::forward_iterator_tag, T >
	{
		/// defining the base type of this iterator
		typedef std::iterator< std::forward_iterator_tag, T >		base_type;
		
		/// defining the own type
		typedef iterator											self_type;
		
		/// defining the type of counter
		typedef std::size_t											size_type;
		
	protected:
		/// value that should stay constant within the iterator whatever happens
		T&															value;
		
		/// variable that count whenever the iterator is incremented
		size_t														counter;
		
	public:
		/** constructor accepting any element of the templated type */
		iterator( T& value_in )
			: base_type()
			, value ( value_in )
			, counter( 0 )
		{ }

		/**
		 * constructor accepting any element of the templated type and an unsigned integer type value.
		 * The unsigned integer values signs the maximum iterations the container should allow.
		 * this can support using the this iterator for the firs two arguments of standard
		 * algorithms like \c std::transform .
		 *
		 * Example use case:\n
		 @code
		  Math::Matrix3x3 valueThatShouldNotChange;
		  Util::identity< Math::Matrix3x3 > idIter( valueThatShouldNotChange, 4 ); // this will stop the algorithm after 4 iterations(increments)
		  std::transform( idIter.begin(), idIter.end(), someOtherInputIterator, AgainSomeOtherOutputIterator, AnyValidFunctionPointer );
		 @endcode
		 *
		 */
		iterator( T& value_in, const std::size_t n_in )
			: base_type()
			, value ( value_in )
			, counter( n_in )
		{ }
		
		/** dereference operator returning reference to the contained object */
		T& operator*()
		{
			return value;
		}
		
		/** pointer operator accessing to access functions of the included object*/
		T* operator->()
		{
			return &value;
		}
		
		self_type& operator+( const self_type& other )
		{
			this->counter += other.counter;
			return *this;
		}
		
		self_type& operator-( const self_type& other )
		{
			this->counter -= other.counter;
			return *this;
		}
		
		/** increment operator for random access operations */
		self_type& operator+( const std::ptrdiff_t inc )
		{
			this->counter += inc;
			return *this;
		}
		
		/** decrement operator for random access operations */
		self_type& operator-( const std::ptrdiff_t dec )
		{
			this->counter -= dec;
			return *this;
		}
		
		self_type& operator+=( const std::ptrdiff_t inc )
		{
			this->counter += inc;
			return *this;
		}
		
		self_type& operator-=( const std::ptrdiff_t dec )
		{
			this->counter -= dec;
			return *this;
		}
		
		/** increment operator, incrementing the internal counter */
		self_type& operator++()
		{
			++counter;
			return *this;
		}
		
		/** decrement operator, decrementing the internal counter */
		self_type& operator--()
		{
			--counter;
			return *this;
		}
		
		/** subscript operator, returning always the contained object no matter which values is requested*/
		template< typename AnyType  >
		T& operator[] ( const AnyType ) const
		{
			return value;
		}
		
		/** inequality operator, used to check if iterator was incremented already enough */
		bool operator!=( const self_type& rhs ) const
		{
			return ( counter != rhs.counter );
		}
		
		/** equality operator, used to check if iterator was incremented already enough */
		bool operator==( const self_type& rhs ) const
		{
			return ( counter == rhs.counter );
		}
		
		/** is lesser operator */
		bool operator<( const self_type& rhs ) const
		{
			return ( counter < rhs.counter );
		}
		
		/** is greater operator */
		bool operator>( const self_type& rhs ) const
		{
			return ( counter > rhs.counter );
		}
	};
	
	typedef const iterator										const_iterator;
	
	/** returning a new iterator pointing always to the contained object*/
	iterator begin()
	{
		return iterator( value, 0 );
	}
	
	/** returning a new iterator pointing always to the contained object */
	iterator end()
	{
		return iterator( value, counter );
	}
	
	/** returning a new iterator pointing always to the contained object*/
	const_iterator cbegin() const
	{
		return const_iterator( value, 0 );
	}
	
	/** returning a new iterator pointing always to the contained object */
	const_iterator cend() const
	{
		return const_iterator( value, counter );
	}
};

} } //namespace Ubitrack::Util

#endif // __UBITRACK_UTIL_IDENTITY_ITERATOR_H__

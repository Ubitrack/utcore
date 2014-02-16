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


// helper structure to have an iterator to always the same element.
// still needs some more functionality like subscript operator, etc.

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
class identity
	: public std::iterator< std::random_access_iterator_tag, T >
{
	typedef identity< T > self_type;
	typedef self_type iterator; 
	
protected:
	const T& value;
	std::size_t n;
	
public:
    identity( const T& value_in )
        : std::iterator< std::random_access_iterator_tag, T >()
		, value ( value_in )
		, n( 0 ){ }

	identity( const T& value_in, const std::size_t n_in )
        : std::iterator< std::random_access_iterator_tag, T >()
		, value ( value_in )
		, n( n_in ){ }
	
	iterator begin()
	{
		return self_type( value, 0 );
	}
	
	iterator end()
	{
		return self_type( value, n );
	}
	
	const T& operator*()
	{
		return value;
	}
	
	const T* operator->()
	{
		return &value;
	}
	
	self_type& operator++()
	{
		++n;
		return *this;
	}
	
	self_type& operator--()
	{
		--n;
		return *this;
	}
	
	T& operator[] ( const std::size_t )
	{
		return value;
	}
	
	bool operator!=( const self_type& rhs ) const
	{
		return ( n != rhs.n );
	}
	
	bool operator==( const self_type& rhs ) const
	{
		return ( n == rhs.n );
	}
};

} } //namespace Ubitrack::Util

#endif // __UBITRACK_UTIL_IDENTITY_ITERATOR_H__

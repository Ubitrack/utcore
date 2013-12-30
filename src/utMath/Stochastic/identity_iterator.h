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

 
 
#ifndef __H__IDENTITY_ITERATOR__
#define __H__IDENTITY_ITERATOR__

#include <utUtil/StaticAssert.h>

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
	typedef typename identity< T > self_type;
	typedef typename self_type iterator; 
	
protected:
	const T value;
	
public:
    identity( const T& value_in )
        : std::iterator< std::random_access_iterator_tag, T >()
		, value ( value_in ){ }
	
	iterator begin()
    {
		return self_type( value );
	}
	
	iterator end()
    {
		UBITRACK_STATIC_ASSERT( (false), THIS_CONTAINER_DOES_NOT_SUPPORT_AN_END_ITERATOR_FUNCTION_ON_PURPOSE );
		return self_type( value );
	}
	
	const T &operator*()
	{
		return (value);
	}
	
	self_type &operator++()
	{
		return *this;
	}
};

} } //namespace Ubitrack::Util

#endif // __H__IDENTITY_ITERATOR__
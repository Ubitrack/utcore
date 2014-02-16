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
 * indices (iterator) template.
 * 
 * This container/iterator adds the possibility to easily
 * use stl-algorithms, like \c std::transform, with a
 * with indexed elements.
 *
 * @todo File will be moved later to another place. -> do not include more
 * often than necessary.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __UBITRACK_UTIL_INDEX_ITERATOR_H__
#define __UBITRACK_UTIL_INDEX_ITERATOR_H__

#include <iterator> // std::iterator, std::input_iterator_tag

namespace Ubitrack { namespace Util { 


// helper structure to have an iterator to always the same element.
// still needs some more functionality like subscript operator, etc.

/**
 * @ingroup util
 * index (iterator) template.
 * 
 * This container/iterator adds the possibility to easily
 * use stl-algorithms, like \c std::transform, with values
 * in a container that are indexed in a seperate container.
 *
 * @tparam value_iterator the type of the iterator pointing to the elements of the container.
 * @tparam index_iterator the type of the iterator pointing to the indices in another container.
 */ 
template< typename value_iterator, typename index_iterator >
class indexed_iterator
	: public std::iterator< std::input_iterator_tag, typename std::iterator_traits< value_iterator >::value_type >
{
	// the base type of this iterator
	typedef typename std::iterator< std::input_iterator_tag, typename std::iterator_traits< value_iterator >::value_type > base_type;
	
	// the type of this iterator itself
	typedef indexed_iterator< value_iterator, index_iterator > self_type;
	
	// the type of the iterator (again)
	typedef self_type iterator;
	
	// the return type (==pointer type) of this iterator
	typedef typename base_type::value_type value_type;
	
	// the type of the indices
	typedef typename std::iterator_traits< index_iterator >::value_type index_type;
	
protected:
	value_iterator itValuesBegin;
	const value_iterator itValuesEnd;
	
	index_iterator pIndex;
	const index_iterator pIndexEnd;
	
	// value of the indices that should be extracted
	const index_type comp_value;
	
public:
	/** Constuctor - expects iterators to the elements and indices */
	indexed_iterator( value_iterator valuesBegin, const value_iterator valuesEnd, index_iterator iBegin, const index_iterator iEnd, const index_type comp_value_in )
        : base_type()
		, itValuesBegin ( valuesBegin )
		, itValuesEnd ( valuesEnd )
		, pIndex( iBegin )
		, pIndexEnd( iEnd )
		, comp_value( comp_value_in )
	{ }
	
	/** Constructor - expects a container to the values and a container to the indices */
	template< typename value_container_type, typename index_container_type >
	indexed_iterator( value_container_type& values, value_container_type& indices, const typename index_container_type::value_type comp_value_in )
        : base_type()
		, itValuesBegin ( values.begin() )
		, itValuesEnd ( values.end() )
		, pIndex( indices.begin() )
		, pIndexEnd( indices.end() )
		, comp_value( comp_value_in )
	{ }
	
	/**  function that estimates the (expected) number of elements */
	const std::size_t size() const
	{
		return std::count( pIndex, pIndexEnd, comp_value );
	}
	
	iterator begin()
	{
		return self_type( itValuesBegin, itValuesEnd, pIndex, pIndexEnd, comp_value );
	}
	
	iterator end()
	{
		return self_type( itValuesBegin, itValuesEnd, pIndex, pIndexEnd, comp_value );
	}
	
	const value_type& operator*()
	{
		return *itValuesBegin;
	}
	
	const value_type* operator->()
	{
		return itValuesBegin;
	}
	
	self_type& operator++()
	{
		const index_iterator it = std::find( pIndex, pIndexEnd, comp_value );
		const std::size_t n = std::distance( pIndex, it );
		std::advance( itValuesBegin, n );
		std::advance( pIndex, n );
		return *this;
	}
	
	// needed to see if iterator reached end
	bool operator!=( const self_type& rhs ) const
	{
		return ( pIndex != rhs.pIndexEnd );
	}
	
	// needed to see if iterator reached end
	bool operator==( const self_type& rhs ) const
	{
		return ( pIndex == rhs.pIndexEnd );
	}
};

} } //namespace Ubitrack::Util

// // warning: do not try this at home! Adding function to std-namesapce :(
// namespace std {

	// template< typename value_iterator, typename index_iterator >
	// distance< Ubitrack::Util::indexed_iterator< value_iterator, index_iterator > >( Ubitrack::Util::indexed_iterator< value_iterator, index_iterator > it1, Ubitrack::Util::indexed_iterator< value_iterator, index_iterator > it2 )
	// {
		// return 0;//it1.size();
	// }
// }

#endif // __UBITRACK_UTIL_INDEX_ITERATOR_H__

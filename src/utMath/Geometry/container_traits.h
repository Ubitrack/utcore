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
 * container traits templates, adding missing functionally to 
 * safe type checking at compilation time.
 * @todo File will be moved later to another place. -> do not include more
 * often than necessary.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

 
 
#ifndef __H__CONTAINER_TRAITS__
#define __H__CONTAINER_TRAITS__

#include <iterator> // iterator_traits

namespace Ubitrack { namespace Util { 


// might be moved later to another location, meant only for internal use
template< typename typeA, typename typeB >
struct is_same
{
	static const bool value = false;
};

// specialization of previous template
template < typename typeA >
struct is_same< typeA, typeA >
{
	static const bool value = true;
};


// Since output iterators (e.g. std::back_insert_itertator via std::back_inserter)
// can be used as well a simple value_type is not enough
// -> therefore designed an own container_traits struct
// might be moved later to another location, meant only for internal use
template< typename iter_type, typename category = std::iterator_traits< iter_type >::iterator_category >
struct container_traits
{
	typedef typename std::iterator_traits< iter_type >::value_type value_type;
};
// specialization of previous template
template< typename iter_type >
struct container_traits< iter_type, std::output_iterator_tag >
{
	typedef typename iter_type::container_type::value_type value_type;
};

}} // namespace Ubitrack::Util

#endif // __H__CONTAINER_TRAITS__
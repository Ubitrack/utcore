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
 * Casts for matrices and vectors.
 * @author Daniel Pustka <pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_UTIL_CAST_ASSIGN_H__
#define __UBITRACK_MATH_UTIL_CAST_ASSIGN_H__

// This include is needed as workaround on Linux, at least with gcc 4.1)
#include <boost/numeric/ublas/functional.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <boost/numeric/ublas/detail/vector_assign.hpp>

namespace Ubitrack { namespace Math { namespace Util {

/**
 * internal class used by the assignments
 * \internal
 */
template< class T1, class T2 >
struct scalar_cast_assign
{
	typedef typename boost::numeric::ublas::type_traits< typename boost::remove_reference< T1 >::type >::value_type cast_type;
	typedef typename boost::numeric::ublas::type_traits< typename boost::remove_reference< T1 >::type >::reference argument1_type;
	typedef typename boost::numeric::ublas::type_traits< T2 >::const_reference argument2_type;
	static const bool computed = false;

	static BOOST_UBLAS_INLINE
	void apply( argument1_type t1, argument2_type t2 ) 
	{ t1 = static_cast< cast_type >( t2 ); }
};


/** 
 * Assigns the contents of m2 to m1, performing an implicit \c static_cast.
 */
template < class T1, class T2 >
void matrix_cast_assign( T1& m1, const T2& m2 )
{ boost::numeric::ublas::matrix_assign< scalar_cast_assign >( m1, m2 ); }


/** 
 * Assigns the contents of v2 to v1, performing an implicit \c static_cast.
 */
template < class T1, class T2 >
void vector_cast_assign( T1& v1, const T2& v2 )
{ boost::numeric::ublas::vector_assign< scalar_cast_assign >( v1, v2 ); }

} } } // namespace Ubitrack::Math::Util

#endif // __UBITRACK_MATH_UTIL_CAST_ASSIGN_H__

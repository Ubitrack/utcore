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
 * Wrapper for simple types.
 * @author Florian Echtler <echtler@in.tum.de>
 */


#ifndef __Scalar_h_INCLUDED__
#define __Scalar_h_INCLUDED__

#include <utCore.h>
#include <boost/serialization/access.hpp>

namespace Ubitrack { namespace Math {


/**
 * @ingroup math
 * Wraps a scalar type so that a Measurement can be derived.
 * This seems overkill, but as this is a simple struct and all
 * member functions are inlined, there's no runtime overhead.
 * @param Builtin a builtin scalar type
 */

template < typename Builtin > 
class Scalar
{
public:

	typedef Builtin value_type;
	
	/** the value of the built-in type*/
	Builtin m_value;
	
	/** Sometimes we need a default constructor */
	inline Scalar()
		: m_value( Builtin() )
	{}

	/**
	 * construct from a built-in type
	 * @param value a scalar
	 */
	inline Scalar( Builtin value )
		: m_value( value )
	{}

	// assignment and copy constructors are automatically generated by the compiler

	/**
	 * conversion operator to cast it back to the built-in type
	 */
	inline operator Builtin() const
	{ return m_value; } 

	
	
protected:
	friend class ::boost::serialization::access;
	
	/** boost::serialization helper */
	template< class Archive >
	void serialize( Archive& ar, const unsigned int )
	{ ar & m_value; }
};

/** comparison for Scalars */
template< class Builtin >
bool operator<( Scalar< Builtin > a, Scalar< Builtin > b )
{ return a.m_value < b.m_value; }

} } // namespace Ubitrack::Math

#endif


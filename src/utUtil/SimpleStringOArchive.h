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
 * @file
 * Conversion of very simple structures to strings.
 * Provides an archive class for boost::serialization
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_UTIL_SIMPLESTRINGOARCHIVE_H_INCLUDED__
#define __UBITRACK_UTIL_SIMPLESTRINGOARCHIVE_H_INCLUDED__

#include <string>
#include <boost/version.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/mpl/bool.hpp>
#include <sstream>
#include <utUtil/Exception.h>

namespace Ubitrack { namespace Util {

/**
 * Writing very simple structures to strings.
 * Provides an archive class for boost::serialization.
 */
class SimpleStringOArchive
{
public:
	/** default constructor */
	SimpleStringOArchive()
		: m_count( 0 )
	{}

	/** get string */
	std::string str()
	{ return m_stream.str(); }
	
protected:
	/** returns the stream used to write to */
	std::ostream& stream()
	{ return m_stream; }
	
	/** pre-write operations */
	void pre()
	{
		if ( m_count++ )
			stream() << ' ';
	}
	
	/** post-write operations */
	void post()
	{}

	/** the string stream to write to */
	std::ostringstream m_stream;

	/** numner of things written */
	unsigned m_count;
	
public:

	/// forward << to &
	template< class T >
	SimpleStringOArchive& operator<<( const T& v )
	{ return *this & v; }

	/// write operator for doubles
	SimpleStringOArchive& operator&( const double& v )
	{ pre(); stream() << v; post(); return *this; }

	/// write operator for floats
	SimpleStringOArchive& operator&( const float& v )
	{ pre(); stream() << v; post(); return *this; }

	/// write operator for ints
	SimpleStringOArchive& operator&( const int& v )
	{ pre(); stream() << v; post(); return *this; }

	/// write operator for unsigned ints
	SimpleStringOArchive& operator&( const unsigned int& v )
	{ pre(); stream() << v; post(); return *this; }

	/// write operator for chars
	SimpleStringOArchive& operator&( const char v )
	{ pre(); stream() << v; post(); return *this; }

	/// write operator for unsigned long longs
	SimpleStringOArchive& operator&( const unsigned long long& v )
	{ pre(); stream() << v; post(); return *this; }

	#if BOOST_VERSION >= 103500
		/// write operator for collection_size_type
		SimpleStringOArchive& operator&( const boost::serialization::collection_size_type& v )
		{ pre(); stream() << v; post(); return *this; }
	#endif
	#if BOOST_VERSION >= 104400
		/// write operator for item_version_type
		SimpleStringOArchive& operator&( const boost::serialization::item_version_type& v )
		{ pre(); stream() << v; post(); return *this; }
	#endif

	/// write operator for all other classes (calls serialize)
	template< class T >
	SimpleStringOArchive& operator&( const T& v )
	{ boost::serialization::serialize( *this, const_cast< T& >( v ), 0 ); return *this; }

	/**
	* ignore binary data
	* bad as a string anyway
	*/
	void save_binary(const void*, size_t) {};
	void load_binary(void*, size_t) {};

	// required for boost::serialization
    typedef boost::mpl::bool_<false> is_loading;
    typedef boost::mpl::bool_<true> is_saving;
	unsigned int get_library_version() const 
	{ return 0; }	
};

} } // namespace Ubitrack::Util

#endif

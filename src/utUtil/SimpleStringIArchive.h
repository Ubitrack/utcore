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
 * Conversion of very simple structures from strings.
 * Provides an archive class for boost::serialization
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_UTIL_SIMPLESTRINGIARCHIVE_H_INCLUDED__
#define __UBITRACK_UTIL_SIMPLESTRINGIARCHIVE_H_INCLUDED__

#include <string>
#include <sstream>
#include <boost/serialization/serialization.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/version.hpp>
#include <utUtil/Exception.h>

namespace Ubitrack { namespace Util {

/**
 * Reading very simple structures from strings.
 * Provides an archive class for boost::serialization.
 */
class SimpleStringIArchive
{
public:
	/** construct archive from string */
	SimpleStringIArchive( const std::string& s )
		: m_string( s )
		, m_stream( m_string )
	{}

protected:
	/** returns the stream used to read from */
	std::istream& stream()
	{ return m_stream; }
	
	/** pre-read operations */
	void pre()
	{
		if ( !m_stream.good() )
			UBITRACK_THROW( "Stream read failure" );
	}
	
	/** post-read operations */
	void post()
	{
		if ( m_stream.fail() )
			UBITRACK_THROW( "Stream read failure" );
	}

	/** the string to store the content */
	std::string m_string;
	
	/** the string stream to read from */
	std::istringstream m_stream;

public:

	/// forward >> to &
	template< class T >
	SimpleStringIArchive& operator>>( T& v )
	{ return *this & v; }

	/// read operator for doubles
	SimpleStringIArchive& operator&( double& v )
	{ pre(); stream() >> v; post(); return *this; }

	/// read operator for floats
	SimpleStringIArchive& operator&( float& v )
	{ pre(); stream() >> v; post(); return *this; }

	/// read operator for ints
	SimpleStringIArchive& operator&( int& v )
	{ pre(); stream() >> v; post(); return *this; }

	/// read operator for unsigned ints
	SimpleStringIArchive& operator&( unsigned int& v )
	{ pre(); stream() >> v; post(); return *this; }

	/// read operator for chars
	SimpleStringIArchive& operator&( char& v )
	{ pre(); stream() >> v; post(); return *this; }

	/// read operator for unsigned long longs
	SimpleStringIArchive& operator&( unsigned long long& v )
	{ pre(); stream() >> v; post(); return *this; }

	#if BOOST_VERSION >= 103500
		/// read operator for collection_size_type
		SimpleStringIArchive& operator&( boost::serialization::collection_size_type& v )
		{ pre(); stream() >> v; post(); return *this; }
	#endif
	#if BOOST_VERSION >= 104400
		/// read operator for item_version_type
		SimpleStringIArchive& operator&( boost::serialization::item_version_type& v )
		{ pre(); stream() >> v; post(); return *this; }
	#endif

	/// let const-nvps through
	template< class T >
	SimpleStringIArchive& operator&( const boost::serialization::nvp< T >& v )
	{ boost::serialization::serialize( *this, const_cast< boost::serialization::nvp< T >& >( v ), 0 ); return *this; }
	
	/// read operator for all other classes (calls serialize)
	template< class T >
	SimpleStringIArchive& operator&( T& v )
	{ boost::serialization::serialize( *this, v, 0 ); return *this; }

	

	/**
	* ignore binary data
	* bad as a string anyway
	*/
	void save_binary(const void*, size_t) {};
	void load_binary(void*, size_t) {};


	// required for boost::serialization
    typedef boost::mpl::bool_<true> is_loading;
    typedef boost::mpl::bool_<false> is_saving;
	unsigned int get_library_version() const 
	{ return 0; }
	void reset_object_address( void*, void* )
	{}	
};

} } // namespace Ubitrack::Util

#endif

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

#include "GlobFiles.h"

// Boost
#include <boost/regex.hpp>

// Ubitrack
#include <utUtil/Exception.h>

#include <log4cpp/Category.hh>
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Util.GlobFiles" ) );

namespace Ubitrack { namespace Util {



void globFiles( const std::string& directory, const std::string & patternString, std::list< boost::filesystem::path >& files, bool globDirectories)
{
	boost::filesystem::path testPath( directory );
	if ( boost::filesystem::is_directory( testPath ) && boost::filesystem::exists( testPath ) )
	{
		boost::regex ext( patternString.c_str() );

		// iterate directory
		boost::filesystem::directory_iterator dirEnd;
		for ( boost::filesystem::directory_iterator it( testPath ); it != dirEnd; it++ )
		{
#ifdef BOOST_FILESYSTEM_I18N
			boost::filesystem::path p( it->path() );
#else
			boost::filesystem::path p( *it );
#endif
			// check for files with suitable extension
			if ( boost::filesystem::exists( p ) && ! boost::filesystem::is_directory( p ) && boost::regex_match( p.filename().string(), ext ) )
			{
				LOG4CPP_TRACE( logger, "Adding file " << p << " to list" );

				files.push_back( p );
			}
			// glob directories, if desired
			else if ( globDirectories && boost::filesystem::exists( p ) && boost::filesystem::is_directory( p ) )
			{
				files.push_back( p );
			}
		}

		// sort, since directory iteration is not ordered on some file systems
		files.sort();

		LOG4CPP_DEBUG( logger, "Sorted list of files" );		
		for ( std::list< boost::filesystem::path >::iterator iter = files.begin(); iter != files.end(); iter ++ )
		{
			LOG4CPP_DEBUG( logger, *iter );		
		}
	}
	else if ( boost::filesystem::exists( testPath ) ) {
		files.push_back( testPath );
	}
	else {
		UBITRACK_THROW( "Invalid path specified" );
	}

	if ( files.size() == 0 ) {
		UBITRACK_THROW( "No suitable files found at the specified location" );	
	}
}

void globFiles( const std::string& directory, const enum FilePattern pattern, std::list< boost::filesystem::path >& files )
{
	std::string patternString = "";
	if ( pattern == PATTERN_DIRECTORIES )
		globFiles ( directory, patternString, files, true);
	else {
		if ( pattern == PATTERN_OPENCV_IMAGE_FILES )
			patternString = ".*\\.(jpg|JPG|png|PNG|bmp|BMP)";
		else if ( pattern == PATTERN_UBITRACK_CALIBRATION_FILES )
			patternString = ".*\\.(cal)";
		else if (pattern == PATTERN_UBITRACK_BOOST_BINARY)
			patternString = ".*\\.(BoostBinary)";
		globFiles ( directory, patternString, files, false);
	}
}

} } // namespace Ubitrack::Util

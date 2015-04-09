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
 * Provides functions to glob files files based on a common path and filename pattern.
 *
 * @author Peter Keitler <peter.keitler@extend3d.de>
 */
 
#ifndef __UBITRACK_UTIL_GLOBFILES_H_INCLUDED__
#define __UBITRACK_UTIL_GLOBFILES_H_INCLUDED__

#include <list>
#include <string>


// Needed until boost version 1.45 because otherwise the deprecated version 2 would be used 
// (see also config/boost where version 2 is currently enforced)
#undef  BOOST_FILESYSTEM_VERSION
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>

#include <utCore.h>	// UBITRACK_EXPORT

namespace Ubitrack { namespace Util {


enum FilePattern {
	PATTERN_OPENCV_IMAGE_FILES,
	PATTERN_UBITRACK_CALIBRATION_FILES,
	PATTERN_DIRECTORIES
};


/**
 * Retrieves all files in the specified directory adhering to the given file name pattern and returns them in the specified list
 */
UBITRACK_EXPORT void globFiles( const std::string& directory, const std::string & patternString, std::list< boost::filesystem::path >& files, bool globDirectories = false);
UBITRACK_EXPORT void globFiles( const std::string& directory, const enum FilePattern pattern, std::list< boost::filesystem::path >& files );


} } // namespace Ubitrack::Util

#endif

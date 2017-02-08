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
 * Provides functions to read and write calibration files.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_UTIL_CALIBFILE_H_INCLUDED__
#define __UBITRACK_UTIL_CALIBFILE_H_INCLUDED__

#include <utMeasurement/Measurement.h> //includes already SharedPtr
#include <utUtil/Exception.h>

#include <string>
#include <fstream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

// get a logger
#include <log4cpp/Category.hh>
static log4cpp::Category& calibLogger( log4cpp::Category::getInstance( "Ubitrack.Utils.CalibFile" ) );

namespace Ubitrack { namespace Util {


/** 
 * read data from a calibration file.
 * @param sFile filename of file to open
 * @param result where the result will be stored
 */
template< typename T >
void readCalibFile( const std::string& sFile, T& result )
{
	// create ifstream
	std::ifstream stream( sFile.c_str() );
	if ( !stream.good() )
		UBITRACK_THROW( "Could not open file " + sFile + " for reading" );

    try {
	    // create iarchive
	    boost::archive::text_iarchive archive( stream );
    	// read data
	    archive >> result;
    } catch (std::exception& ) {
        UBITRACK_THROW( "Could not read ubitrack file" );
    }
}

/**
* read data from a binary calibration file.
* @param sFile filename of file to open
* @param result where the result will be stored
*/
template< typename T >
void readBinaryCalibFile(const std::string& sFile, T& result)
{
	// create ifstream
	std::ifstream stream(sFile.c_str(), std::ios::in | std::ios::binary);
	if (!stream.good())
		UBITRACK_THROW("Could not open file " + sFile + " for reading");

	try {
		// create iarchive
		boost::archive::binary_iarchive archive(stream);
		// read data
		archive >> result;
	}
	catch (std::exception& e) {
		std::string tmp = e.what();
		UBITRACK_THROW("Could not read ubitrack file: " + tmp);
	}
}


/** 
 * read data from a calibration file (specialized for Measurements).
 * @param sFile filename of file to open
 * @param result where the result will be stored
 */
template< typename T >
void readCalibFile( const std::string& sFile, Measurement::Measurement< T >& result )
{
	// initialize measurement
	result = Measurement::Measurement< T >( 0, boost::shared_ptr< T >( new T() ) );

	// create ifstream
	std::ifstream stream( sFile.c_str() );
	if ( !stream.good() )
		UBITRACK_THROW( "Could not open file " + sFile + " for reading" );

	// create iarchive
	boost::archive::text_iarchive archive( stream );

    try {
	    // read data
	    archive >> result;
    } catch (std::exception& ) {
        UBITRACK_THROW( "Wrong file format" );
    }
}

/** 
 * Read data from a calibration file (containing Measurements)
 * but drops the timestamp and just initializes the payload.
 * This is a pure temporary helper function until no calibration
 * files has an timestamp anymore.
 * @param sFile filename of file to open
 * @param result Where the result will be stored
 */
template< typename T >
void readCalibFileDropMeasurement( const std::string& sFile, T& result )
{
    LOG4CPP_WARN( calibLogger, "Reading calibration files with measurement overhead. Consider using files without timestamp!" );
	// initialize measurement
	Measurement::Measurement< T > interResult ( 0, boost::shared_ptr< T >( new T() ) );

	// create ifstream
	std::ifstream stream( sFile.c_str() );
	if ( !stream.good() )
		UBITRACK_THROW( "Could not open file " + sFile + " for reading" );

	// create iarchive
	boost::archive::text_iarchive archive( stream );

	try {
	    // read data
	    archive >> interResult;
    } catch (std::exception& ) {
        UBITRACK_THROW( "Wrong file format" );
    }

    // drop measurement and write to out param
    result = *interResult;
}

/** 
 * write data to a calibration file.
 * @param sFile filename of file to open
 * @param data data to save
 */
template< typename T >
void writeCalibFile( const std::string& sFile, const T& data )
{
	// create ofstream
	std::ofstream stream( sFile.c_str() );
	if ( !stream.good() )
		UBITRACK_THROW( "Could not open file " + sFile + " for writing" );

	// create iarchive
	boost::archive::text_oarchive archive( stream );

	// read data
	archive << data;
}

/**
* write data to a calibration file.
* @param sFile filename of file to open
* @param data data to save
*/
template< typename T >
void writeBinaryCalibFile(const std::string& sFile, const T& data)
{
	// create ofstream
	std::ofstream stream(sFile.c_str(), std::ios::out | std::ios::binary);
	if (!stream.good())
		UBITRACK_THROW("Could not open file " + sFile + " for writing");

	// create iarchive
	boost::archive::binary_oarchive archive(stream);

	// read data
	archive << data;
}
} } // namespace Ubitrack::Util

#endif


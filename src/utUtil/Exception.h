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
 * Contains standard runtime exceptions for Ubitrack library.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_UTIL_EXCEPTION_H_INCLUDED__
#define __UBITRACK_UTIL_EXCEPTION_H_INCLUDED__

#include <stdexcept>
#include <string>
#include <iostream>
#include <utCore.h>

namespace Ubitrack { namespace Util {

/**
 * Base class for all ubitrack exceptions.
 * Adds file and line number information to the c++ standard exceptions.
 */
class UBITRACK_EXPORT Exception 
	: public std::runtime_error
{
public:
	/** 
	 * Constructs the message. 
	 *
	 * @param sMessage message to show to the user
	 * @param nLine line number where the exception was thrown.
	 * @param sFile source file where the exception was thrown.
	 */
	Exception( const std::string& sMessage, unsigned nLine = 0, const char* sFile = NULL );

	/** the default destructor does not specify throw(), so we have to do it by hand */
	~Exception() throw();
	
	/** returns the line number where the exception was thrown. */
	unsigned line() const
	{ return m_nLine; }
	
	/** returns the source file where the exception was thrown. */
	const std::string& file() const
	{ return m_sFile; }
	
protected:
	unsigned m_nLine;
	std::string m_sFile;
};


/** Streaming operator for ubitrack exceptions. */
UBITRACK_EXPORT std::ostream& operator<<( std::ostream& o, const Exception& e );


} } // namespace Ubitrack::Util


/** 
 * Macro to throw an exception.
 * File and line number will we generated automatically.
 */
#define UBITRACK_THROW( message ) throw Ubitrack::Util::Exception( message, __LINE__, __FILE__ )

#endif

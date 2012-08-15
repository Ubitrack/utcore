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
 * Implements standard runtime exceptions for Ubitrack library.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 
 
#include "Exception.h"
#include <log4cpp/Category.hh>
 
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Util.Exception" ) );

namespace Ubitrack { namespace Util {


Exception::Exception( const std::string& sMessage, unsigned nLine, const char* sFile )
	: std::runtime_error( sMessage )
	, m_nLine( nLine )
	, m_sFile( sFile ? sFile : "" )
{
	LOG4CPP_DEBUG( logger, "Exception thrown in " << sFile << ":" << nLine << ", message: " << sMessage );
}

/* bug fix/workaround for a problem with VS2005 which has problems when the exception class is
 * not UBITRACK_EXPORT and for this, some method must be in an implementation file.
 */
Exception::~Exception() throw()
{}


std::ostream& operator<<( std::ostream& o, const Exception& e )
{
	LOG4CPP_TRACE( logger, "Exception \"" << e.what() << "\" from " << e.file() << ":" << e.line() );

	o << "Exception \"" << e.what() << "\" from " << e.file() << ":" << e.line();
	return o;
}

} } // namespace Ubitrack::Util

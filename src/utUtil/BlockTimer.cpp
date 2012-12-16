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
 * Implementation of a high resolution timer to measure execution time of code blocks
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include <iostream>
#include <sstream>
#include <iomanip>

#include "BlockTimer.h"

namespace Ubitrack { namespace Util { 


BlockTimer::~BlockTimer()
{
	if ( m_pLogger && m_nRuns ) 
	{ 
		std::ostringstream s;
		s << *this;
		m_pLogger->log( log4cpp::Priority::debug, s.str(), m_sCodeFile.c_str(), m_nCodeLine );
	}
}


void BlockTimer::initializeStart( const char* sCodeFile, unsigned nCodeLine )
{ 
	m_sCodeFile = sCodeFile; 
	m_nCodeLine = nCodeLine;
}


void BlockTimer::initializeEnd()
{
	m_bInitialized = true;
	// TODO: add automatic hierarchy detection
}




std::ostream& operator<<( std::ostream& s, const BlockTimer& t )
{		
	unsigned long long totalRunTime = getHighPerformanceCounter() - t.m_startTime;
	return s << std::setw( 30 ) << t.getName()
		<< " runs: " << std::setw( 6 ) << t.getRuns()
		<< ", total: " << std::setw( 7 ) << t.getTotalTime()
		<< "ms, avg: " << std::setw( 7 ) << t.getAvgTime() << "ms"
		<< "ms, totalRuntime: " << std::setw( 7 ) << totalRunTime << "ms"
		<< "ms, call per second: " << std::setw( 7 ) << totalRunTime/1000.0 /  t.getRuns() << "ms";
		
}

} } // namespace Ubitrack::Util

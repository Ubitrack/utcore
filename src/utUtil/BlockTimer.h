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
 * Header file for a high resolution timer to measure execution time of code blocks
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#include <string>
#include <iostream>
#include <log4cpp/Category.hh>
#include <boost/utility.hpp>
#include <utCore.h>
#include <utUtil/OS.h>

namespace Ubitrack { namespace Util {

// forward decl
class BlockTimer;


/** stream output operator for BlockTimers */
UBITRACK_EXPORT std::ostream& operator<<( std::ostream& s, const BlockTimer& t );


/**
 * Times a block of execution.
 *
 * The timer is started by instantiating a \c Time object and stopped when it leaves scope. 
 * The result of multiple runs is summed up. The BlockTimer result can be directly printed to an ostream.
 */
class UBITRACK_EXPORT BlockTimer
{
public:

	/** times a block of execution */
	class Time
		: public boost::noncopyable
	{
	public:
		/** 
		 * Starts timing of a block of execution.
		 * @param rTiming BlockTimer object in which to store the result
		 */
		Time( BlockTimer& rTimer )
			: m_rTimer( rTimer )
			, m_startTime( getHighPerformanceCounter() )
		{}
		
		/** 
		 * Starts timing of a block of execution and sets the code location if not already set.
		 * Internally used by the UBITRACK_TIME macro.
		 * @param rTiming BlockTimer object in which to store the result
		 * @param sCodeFile code location filename 
		 * @param nCodeLine code location line
		 */
		Time( BlockTimer& rTimer, const char* sCodeFile, unsigned nCodeLine )
			: m_rTimer( rTimer )
			, m_startTime( getHighPerformanceCounter() )
		{
			if ( !m_rTimer.initialized() )
				m_rTimer.initializeStart( sCodeFile, nCodeLine );
		}
		
		/**
		 * Stops execution and stores the result in the BlockTimer object given in the constructor
		 */
		~Time()
		{ 
			m_rTimer.addMeasurement( getHighPerformanceCounter() - m_startTime ); 
			
			if ( !m_rTimer.initialized() )
				m_rTimer.initializeEnd();
		}

	protected:
		BlockTimer& m_rTimer;
		unsigned long long m_startTime;
	};

	/** 
	 * constructs and empty block timer object
	 * @param sName name of the timer (for display)
	 * @param sLoggingCategory log4cpp category to which to print the result when the timer object leaves scope
	 */
	BlockTimer( const std::string& sName, const std::string& sLoggingCategory = std::string() )
		: m_sName( sName )
		, m_pLogger( sLoggingCategory.empty() ? 0 : &log4cpp::Category::getInstance( sLoggingCategory ) )
		, m_bInitialized( false )
		, m_nRuns( 0 )
		, m_nTicks( 0 )
	{}

	/** 
	 * constructs and empty block timer object
	 * @param sName name of the timer (for display)
	 * @param logger log4cpp logger to which to print the result when the timer object leaves scope
	 */
	BlockTimer( const std::string& sName, log4cpp::Category& logger )
		: m_sName( sName )
		, m_pLogger( &logger )
		, m_bInitialized( false )
		, m_nRuns( 0 )
		, m_nTicks( 0 )
	{}

	/** destructor, prints result if a stream was given to the constructor */
	~BlockTimer();
	
	/** adds a timer run to the internal state */
	void addMeasurement( unsigned long long ticks )
	{ m_nRuns++; m_nTicks += ticks; }
	
	const std::string& getName() const
	{ return m_sName; }
	
	/** returns the total time in ms */
	double getTotalTime() const
	{ return m_nTicks / getHighPerformanceFrequency() * 1000; }
	
	/** returns the average time in ms */
	double getAvgTime() const
	{ return m_nTicks / ( getHighPerformanceFrequency() * m_nRuns ) * 1000; }
	
	/** returns the number of times the timer was run */
	unsigned getRuns() const
	{ return m_nRuns; }
	
	/** Are the additional informations about the timer initialized? */
	bool initialized() const
	{ return m_bInitialized; }
	
	/** Initialization at begin of first run. Sets the code location (for more useful output) */
	void initializeStart( const char* sCodeFile, unsigned nCodeLine );
		
	/** Initialization at end of first run */
	void initializeEnd();
		
protected:
	std::string m_sName;
	log4cpp::Category* m_pLogger;
	std::string m_sCodeFile;
	unsigned m_nCodeLine;
	bool m_bInitialized;
	
	unsigned m_nRuns;
	unsigned long long m_nTicks;
};


#ifndef UBITRACK_NOTIME
	
	/**
	 * Macro for convenient timing.
	 * Add to beginning of code block whose execution time you want to measure.
	 * @param timer BlockTimer for aggregation of results
	 */
	#define UBITRACK_TIME( timer ) Ubitrack::Util::BlockTimer::Time _timeVar_##timer( timer, __FILE__, __LINE__ )
#else
	#define UBITRACK_TIME( timer )
#endif

 
 } } // namespace Ubitrack::Util

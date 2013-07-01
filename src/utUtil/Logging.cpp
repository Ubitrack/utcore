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
 * Implements Ubitrack-specific functions related to the logging framework.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable: 4290 )
#endif

#ifdef ANDROID
#include <log4cpp/AndroidLogAppender.hh>
#else 
#include <log4cpp/OstreamAppender.hh>
#endif



#include <log4cpp/PatternLayout.hh>
#include <log4cpp/Category.hh>
#include <log4cpp/PropertyConfigurator.hh>

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include <iostream>
#include "Logging.h"



namespace Ubitrack { namespace Util {

// Initializes the logger
void initLogging( const char* sConfigFile )
{
	#ifdef ANDROID
	log4cpp::Appender* app = new log4cpp::AndroidLogAppender( "stderr");
	log4cpp::PatternLayout* layout = new log4cpp::PatternLayout();
	layout->setConversionPattern( "%d{%H:%M:%S.%l} %6p %20f:%-3l %m   (%c)%n" );
//	layout->setConversionPattern( "%R %p %c %x: %m%n" );
	app->setLayout( layout );

	log4cpp::Category::getRoot().setAdditivity( false );
	log4cpp::Category::getRoot().addAppender( app );
	log4cpp::Category::getRoot().setPriority( log4cpp::Priority::INFO ); // default: INFO
	log4cpp::Category::getInstance( "Ubitrack.Events" ).setPriority( log4cpp::Priority::NOTICE ); // default: NOTICE

	#else
	// try to initialize logging from file
	try
	{ 
		log4cpp::PropertyConfigurator::configure( sConfigFile ); 
		return;
	}
	catch ( ... )
	{ } // continue with the default logging settings

	// add a stderr appender with some nice layout and set priority to INFO and event priority to notice
	log4cpp::Appender* app = new log4cpp::OstreamAppender( "stderr", &std::cerr );
	log4cpp::PatternLayout* layout = new log4cpp::PatternLayout();
	layout->setConversionPattern( "%d{%H:%M:%S.%l} %6p %20f:%-3l %m   (%c)%n" );
//	layout->setConversionPattern( "%R %p %c %x: %m%n" );
	app->setLayout( layout );

	log4cpp::Category::getRoot().setAdditivity( false );
	log4cpp::Category::getRoot().addAppender( app );
	log4cpp::Category::getRoot().setPriority( log4cpp::Priority::INFO ); // default: INFO
	log4cpp::Category::getInstance( "Ubitrack.Events" ).setPriority( log4cpp::Priority::NOTICE ); // default: NOTICE
	#endif
}

} } // namespace Ubitrack::Util

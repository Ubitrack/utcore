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


/* 
   This file is only for doxygen content.
   No Code goes in here.
*/

/*
   Group definitions
*/

/** @defgroup datastructures Datastructures
 *  This group contains datastructrues used by the Ubitrack library.
 */
/** @defgroup tracking_algorithms Tracking and Calibration Algorithms
 *  This group contains general tracking and calibration algorithms
 *  for standalone use and for use in dataflow components.
 */
/** @defgroup math General purpose mathematics functions
 *  This group contains general mathematical implementations which
 *  can be used in a wider range of applications.
 */


/*
   Main page content
*/

/**
 * @mainpage Ubitrack Library Overview
 * @section Introduction
 * The Ubitrack Library is a lightweight and efficient tracking library. Its focus is
 * on spatial relationship patterns and sensor fusion.
 * @section Usage
 * example.cpp:
 * @code
#include <Ubitrack/Facade/AdvancedFacade.h>

void callback( const Ubitrack::Measurement::Pose& pose )
{
   std::cout << pose << std::endl;
}

int main( int argc, char* argv[] )
{
  // initialize ubitrack library
  Ubitrack::Facade::AdvancedFacade utFacade;
  utFacade.loadDataflow( "example.utql" );

  // set callback(s) on ApplicationPushSink
  utFacade.setCallback< Ubitrack::Measurement::Pose >( "Sink1", &callback );

   // start tracking
  utFacade.startDataflow();

  // Now do the application-specific things
  // [...]

  // stop tracking
  utFacade.stopDataflow();

  return 0;
}
@endcode
 * example.utql
 * @code
<!--?xml version="1.0" encoding="utf-8"?-->
<UTQLResponse>

<Pattern name="Art6D" id="Art1">
	<Output>
		<Node name="Art" id="Art">
			<Attribute name="artPort" value="5000"/>
		</Node>
		<Node name="Body" id="Body1"/>
		<Edge name="ArtToTarget" source="Art" destination="Body">
			<Attribute name="artBodyId" value="3"/>
			<Attribute name="artType" value="6d"/>
			<Attribute name="type" value="6D"/>
			<Attribute name="mode" value="push"/>
		</Edge>
	</Output>

	<DataflowConfiguration>
		<UbitrackLib class="ArtTracker"/>
	</DataflowConfiguration>
</Pattern>

<Pattern name="ApplicationPushSinkPose" id="Sink1">
	<Input>
		<Node name="A" id="Art"/>
		<Node name="B" id="Body1"/>
		<Edge name="Input" source="A" destination="B" pattern-ref="Art1" edge-ref="ArtToTarget"/>
	</Input>

	<DataflowConfiguration>
		<UbitrackLib class="ApplicationPushSinkPose"/>
	</DataflowConfiguration>
</Pattern>

</UTQLResponse>
 * @endcode
 */
 
/*
   Namespace documentation
*/

/**
 * High-level functions for end-user applications
 */

namespace Ubitrack {

	/**
	 * Calibration algorithms
	 */

	namespace Calibration {}

	/**
	 * Mathematical types and functions
	 */

	namespace Math {}

	/**
	 * Tracking measurement datastructures
	 */

	namespace Measurement {}

	/**
	 * Miscellaneous helper functions
	 */

	namespace Util {}

}


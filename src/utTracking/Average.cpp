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
 * @ingroup dataflow_components
 * @file
 * This class calculates the average of a List measurement.
 *
 * @author Florian Echtler <echtler@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */

// Ubitrack
#include "Average.h"
/*
namespace Ubitrack { namespace Tracking {
	typedef Average< Math::Vector< double, 3 > , Math::Vector< double, 3 >     >  AveragePosition;
	typedef Average< Measurement::Distance, Measurement::Distance > >  DistanceListAverage;
	typedef Average< Measurement::Position2D, Measurement::Position2D > > PositionList2DAverage;
	typedef Average< Measurement::Position, Measurement::Position > > PositionListAverage;
	typedef Average< Measurement::Pose, Measurement::Pose > > PoseListAverage;
	typedef Average< Measurement::Rotation, Measurement::Rotation > > RotationListAverage;
	
	// with Covariance
	typedef Average< Measurement::Position, Measurement::ErrorPosition > > PositionListAverageError;
	typedef Measurement::Pose, Measurement::ErrorPose > > PoseListAverageError;



} } // namespace Ubitrack::Components

*/
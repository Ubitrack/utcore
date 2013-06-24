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
 * @ingroup tracking_algorithms
 * @file
 * Implements functions for bundle adjustment.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#include "BundleAdjustment.h"
#ifdef HAVE_LAPACK

#include "Function/BundleAdjustment.h"

// get a logger
#include <log4cpp/Category.hh>
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Calibration.BundleAdjustment" ) );

//#define OPTIMIZATION_LOGGING
//static log4cpp::Category& optLogger( log4cpp::Category::getInstance( "Ubitrack.Calibration.BundleAdjustment" ) );

#include <utMath/GaussNewton.h>
#include <utMath/LevenbergMarquardt.h>

// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;


namespace Ubitrack { namespace Calibration {

/** 
 * \internal
 * performs the bundle adjustment according to the network description 
 */
template< class T >
T bundleAdjustmentImpl( BundleAdjustmentNetwork< T >& net )
{
	// build problem class
	BundleAdjustmentFunction< T > func( net );

	// build measurement vector
	ublas::vector< T > measurement( func.measurementSize() );
	func.buildMeasurementVector( measurement );

	// build parameter vector
	ublas::vector< T > parameters( func.parameterSize() );
	func.buildParameterVector( parameters );

	LOG4CPP_DEBUG( logger, "Optimizing " << parameters.size() << " parameters using " << measurement.size() << " measurements" );

	// perform optimization
	T residual = Math::levenbergMarquardt( func, parameters, measurement, Math::OptTerminate( 200, 1e-6 ), Math::OptNoNormalize() );

	// update network description
	func.updateParametersFromVector( parameters );

	return residual;
}


/** 
 * performs the bundle adjustment according to the network description 
 */
float bundleAdjustment( BundleAdjustmentNetwork< float >& net )
{ 
	return bundleAdjustmentImpl( net ); 
}

/** 
 * performs the bundle adjustment according to the network description 
 */
double bundleAdjustment( BundleAdjustmentNetwork< double >& net )
{ 
	return bundleAdjustmentImpl( net ); 
}

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

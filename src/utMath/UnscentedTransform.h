/*
 * Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2008, Technische Universitaet Muenchen, and individual
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
 * @author Tobias Reichl <reichl@in.tum.de>
 */

#ifndef __UBITRACK_MATH_UNSCENTEDTRANSFORM_H_INCLUDED__
#define __UBITRACK_MATH_UNSCENTEDTRANSFORM_H_INCLUDED__
 

#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/LevenbergMarquardt.h>

// #include <boost/numeric/ublas/matrix_proxy.hpp>
// #include <boost/numeric/ublas/vector_proxy.hpp>

namespace Ubitrack { namespace Math {

/**
 * Performs an Unscented Transform based on a set of measurements in 2D, a given variance (the probability distribution
 * in 2D is assumed to be isotrophic)  a 2D->6D function, and returns the predicted 6D covariance.
 */

// Note: the UnscentedTransform might be generalizable to arbitrary measurements and matrix dimensions,
//    the covariance computation is specialized here (reduction from seven to six DOF)...

template< class PType, class VType >
Math::Matrix< 6, 6 > unscentedTransform(
		const std::vector< Math::Vector< 2, VType > >& measurements,
		VType variance, PType& problem)
{
	namespace ublas = boost::numeric::ublas;

	// Get standard deviation
	VType stddev = sqrt(variance);
	
	// Storage for optimized parameters
	std::vector< Math::Vector< 0, VType > > optimizedParameters;

	// First guess of parameters
	// TODO: better values than zero?
	Math::Vector< 7, VType > params = Math::Vector< 7, VType >::zeros();
	// Avoid degenerate (absolute value of zero) quaternion
	params( 4 ) = 1;

	// Combine measurements into single vector
	Math::Vector< 0, VType > measurementsCombined ( 2 * measurements.size() );
	for ( unsigned i = 0; i < measurements.size(); i++ )
	{
		measurementsCombined[ 2 * i ] = measurements[ i ][ 0 ];
		measurementsCombined[ 2 * i + 1 ] = measurements[ i ][ 1 ];
	}
	
	// Create undisturbed set
	Math::Vector< 0, VType > sigmaSet ( 2 * measurements.size() );
	sigmaSet = measurementsCombined;
	
	// LevenbergMarquadt on undisturbed set
	VType residual = Math::levenbergMarquardt(problem, params, sigmaSet,
		Math::OptTerminate( 200, 1e-6 ), Math::OptNoNormalize());

	// Store optimized parameters
	optimizedParameters.push_back(params);
	
	// Iterate over all measurements
	for ( unsigned i = 0; i < measurementsCombined.size(); i++ )
	{
		// Copy undisturbed set and disturb it
		sigmaSet = measurementsCombined;
		sigmaSet[ i ] += stddev;

		// Reset parameters
		params = Math::Vector< 7, VType >::zeros();
		params( 4 ) = 1;

		// LevenbergMarquadt on disturbed set
		residual = Math::levenbergMarquardt(problem, params, sigmaSet,
			Math::OptTerminate( 200, 1e-6 ), Math::OptNoNormalize());

		// Debugging: output progress
		// std::cout << "+";

		// Store optimized parameters
		optimizedParameters.push_back(params);

		// Copy undisturbed set and disturb it
		sigmaSet = measurementsCombined;
		sigmaSet[ i ] -= stddev;

		// Reset parameters
		params = Math::Vector< 7, VType >::zeros();
		params( 4 ) = 1;

		// LevenbergMarquadt on disturbed set
		residual = Math::levenbergMarquardt(problem, params, sigmaSet,
			Math::OptTerminate( 200, 1e-6 ), Math::OptNoNormalize());

		// Debugging: output progress
		// std::cout << "-";

		// Store optimized parameters
		optimizedParameters.push_back(params);
	}
	// Debugging: output progress
	// std::cout << std::endl;
	
	// Compute average pose
	Math::Vector< 7, double > avgPose( Math::Vector< 7, double >::zeros(); );
	for ( typename std::vector< Math::Vector< 0, VType > >::iterator it = optimizedParameters.begin();
			it != optimizedParameters.end(); it++ )
	{
		if ( it != optimizedParameters.begin() && 
			ublas::inner_prod( ublas::subrange( *it, 3, 7 ), ublas::subrange( avgPose, 3, 7 ) ) < 0 )
			ublas::subrange( *it, 3, 7 ) *= -1;
		avgPose += *it;
	}
	avgPose /= optimizedParameters.size();
	
	// Debugging: output average pose
	// std::cout << "Average pose:" << std::endl;
	// std::cout << avgPose << std::endl;

	// Get average rotation
	Math::Quaternion avgQuat( Math::Quaternion::fromVector( ublas::subrange( avgPose, 3, 7 ) ) );
	avgQuat.normalize();
	
	// Compute covariance
	Math::Matrix< 6, 6 > covariance( Math::Matrix< 6, 6 >::zeros() );
	for ( typename std::vector< Math::Vector< 0, VType > >::iterator it = optimizedParameters.begin();
			it != optimizedParameters.end(); it++ )
	{
		Math::Vector< 6 > localError;
		ublas::subrange( localError, 0, 3 ) = ublas::subrange( *it, 0, 3 ) - ublas::subrange( avgPose, 0, 3 );;
		
		Math::Quaternion qLocal((*it)[ 3 ], (*it)[ 4 ], (*it)[ 5 ], (*it)[ 6 ]);
		qLocal.normalize();

		Math::Quaternion qDiff( avgQuat * ~qLocal );
		localError( 3 ) = qDiff.x();
		localError( 4 ) = qDiff.y();
		localError( 5 ) = qDiff.z();
		if ( qDiff.w() < 0 )
			ublas::subrange( localError, 3, 6 ) *= -1;
		covariance += ublas::outer_prod( localError, localError );

		// Debugging: output error inner product
		// VType innerProd = ublas::inner_prod( localError, localError );
		// std::cout << "Inner product:" << std::endl;
		// std::cout << innerProd << std::endl;
	}
	covariance /= optimizedParameters.size();
	
	// Return 6D covariance
	return covariance;

}
	
} } // namespace Ubitrack::Math

#endif

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
 * @ingroup math
 * @file
 * Helper classes common to all optimizers
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 
 
#ifndef __UBITRACK_MATH_OPTIMIZATION_H_INCLUDED__
#define __UBITRACK_MATH_OPTIMIZATION_H_INCLUDED__

#include <math.h>

// to turn on logging of internal processing, create a log4cpp::Category object called "optLogger"
// and #define OPTIMIZATION_LOGGING before including this header 
#ifdef OPTIMIZATION_LOGGING
	#include <boost/numeric/ublas/io.hpp>
	#define OPT_LOG_TRACE( message ) LOG4CPP_TRACE( optLogger, message )
	#define OPT_LOG_DEBUG( message ) LOG4CPP_DEBUG( optLogger, message )
	#define OPT_LOG_INFO( message ) LOG4CPP_INFO( optLogger, message )
#else
	#define OPT_LOG_TRACE( message ) 
	#define OPT_LOG_DEBUG( message ) 
	#define OPT_LOG_INFO( message ) 
#endif


namespace Ubitrack { namespace Math { 

/** 
 * Termination criteria that makes the optimizer terminate after n iterations or until the 
 * residual change is smaller than some percentage, whichever comes first. 
 */
class OptTerminate
{
public:
	/**
	 * Constructor.
	 * @param maxIterations maximum number of iterations, <= 0 for unlimited iterations
	 * @param precision stops if the residual r changes from one iteration to the next by less 
	 *    than r*precision.
	 */
	OptTerminate( const std::size_t maxIterations, double precision = 0.0 )
		: m_maxIterations( maxIterations )
		, m_precision( precision )
	{}

	/** this function is evaluated by the optimizer */
	bool operator()( const std::size_t iterations, const double resPrev, const double resNow ) const
	{ 
		return ( m_maxIterations > 0 && iterations >= m_maxIterations ) ||
			( m_precision != 0.0 && fabs( resPrev - resNow ) < m_precision * resNow );
	}

protected:
	const std::size_t m_maxIterations;
	const double m_precision;
};


/** default normalization, which does nothing */
struct OptNoNormalize
{
	template< class X, class Y > 
	void evaluate( X&, const Y& ) const
	{}
};


class OptNoWeightFunction
{
public:
	bool noWeights() const
	{ return true; }

	template< class X, class Y > 
	void computeWeights( const X& errorVector, Y& weightVector ) const
	{}
};

} } // namespace Ubitrack::Math

#endif

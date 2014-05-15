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
 * Functions for ransac tooltip/hotspot calibration.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGROITHM_TOOLTIP_RANSAC_H_INCLUDED__
#define __UBITRACK_ALGROITHM_TOOLTIP_RANSAC_H_INCLUDED__


#include "LeastSquare.h" // includes std::vector/Pose
#include <utCore.h>
#include <utMath/Blas1.h>
#include <utMath/Optimization/Ransac.h>

namespace Ubitrack { namespace Algorithm { namespace ToolTip {

/**
 * @internal function object that provides estimation and evaluation
 * functions for a ransac tooltip calibration.
 */
template< typename T >
struct Ransac
{
public:

	typedef T value_type;
	
	typedef Math::Vector< T, 6 > result_type;
	
	/**
	 * @internal computes a to the tooltip calibration from given parameters
	 */
	struct Estimator
	{
		public:
		
		template< typename InputIterator >
		bool operator()( Math::Vector< T, 6 >& resultVector, const InputIterator iBegin, const InputIterator iEnd ) const
		{
			Math::Vector< T, 3 > pw;
			Math::Vector< T, 3 > pm;
			if ( !estimatePosition3D_6D( pw, iBegin, iEnd, pm ) )
				return false;
				
			resultVector[ 0 ] = pw[ 0 ];
			resultVector[ 1 ] = pw[ 1 ];
			resultVector[ 2 ] = pw[ 2 ];
			resultVector[ 3 ] = pm[ 0 ];
			resultVector[ 4 ] = pm[ 1 ];
			resultVector[ 5 ] = pm[ 2 ];
			return true;
		}
	};
	
	/**
	 * @internal computes euclidean distance of transformed point to original point
	 */
	struct Evaluator
	{
		public:
		
		T operator()( const Math::Vector< T, 6 > &toolTip, const Math::Pose &p ) const
		{
			const Math::Vector< T, 3 > tip = p * Math::Vector< T, 3 > ( toolTip[ 3 ], toolTip[ 4 ], toolTip[ 5 ] );
			const Math::Vector< T, 3 > offset( toolTip[ 0 ], toolTip[ 1 ], toolTip[ 2 ] );
			const Math::Vector< T, 3 > diff( offset - tip );
			return Ubitrack::Math::norm_2( diff );
		}
	};
};

template< typename T, typename InputIterator >
bool estimatePosition3D_6D( Math::Vector< T, 3 >& pw
	, const InputIterator iBegin
	, const InputIterator iEnd
	, Math::Vector< T, 3 >& pm
	, const Math::Optimization::RansacParameter< T >& params )
{
	Math::Vector< T, 6 > resultVector;
	const std::size_t inlier = Math::Optimization::ransac( iBegin, iEnd, resultVector, ToolTip::Ransac< T >(), params  );
	
	pw = Math::Vector< T, 3 > ( resultVector[ 0 ], resultVector[ 1 ], resultVector[ 2 ] );
	pm = Math::Vector< T, 3 > ( resultVector[ 3 ], resultVector[ 4 ], resultVector[ 5 ] );
	
	return ( inlier > 0 );
}

UBITRACK_EXPORT bool estimatePosition3D_6D( Math::Vector3f& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3f& pm
	, const Math::Optimization::RansacParameter< float >& params );
	

UBITRACK_EXPORT bool estimatePosition3D_6D( Math::Vector3d& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector3d& pm
	, const Math::Optimization::RansacParameter< double >& params );	

}}} // namespace Ubitrack::Algorithm::ToolTip

#endif //__UBITRACK_ALGROITHM_TOOLTIP_RANSAC_H_INCLUDED__
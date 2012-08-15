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
 * @ingroup calibration
 * @file
 * functions to normalize a pose computed from projection
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_PROJECTIVEPOSENORMALIZE_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_PROJECTIVEPOSENORMALIZE_H_INCLUDED__
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Function that normalizes the quaternion of a given pose and factors the 
 * difference into the translation. This improves convergence of the 
 * optimization for single-camera pose estimation.
 *
 * The pose is given as a 7-vector (tx, ty, tz, qx, qy, qz, qw)
 */
struct ProjectivePoseNormalize
{
	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 7; }

	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		namespace ublas = boost::numeric::ublas;

		// compute length of quaternion
		typename VT1::value_type fQuatLenSq = ublas::inner_prod( ublas::subrange( input, 3, 7 ), ublas::subrange( input, 3, 7 ) );
		typename VT1::value_type fQuatLen = std::sqrt( fQuatLenSq );

		// normalize quaternion
		ublas::subrange( result, 3, 7 ) = ublas::subrange( input, 3, 7 ) / fQuatLen;
		
		// scale translation
		ublas::subrange( result, 0, 3 ) = ublas::subrange( input, 0, 3 ) / fQuatLenSq;
	}

	// jacobians not needed so far
};

} } } // namespace Ubitrack::Calibration::Function

#endif

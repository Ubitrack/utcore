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
 * Functions for homography estimation.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_CALIBRATION_HOMOGRAPHY_H_INCLUDED__
#define __UBITRACK_CALIBRATION_HOMOGRAPHY_H_INCLUDED__



#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/Matrix.h>
#include <utMath/Vector.h>
#include <vector>

namespace Ubitrack { namespace Calibration {


/**
 * @ingroup tracking_algorithms
 * Computes a general homography using a linear DLT method.
 *
 * The result is a homography H that maps points x to points x' via x'=Hx.
 * See Hartley&Zisserman for details.
 *
 * Note: also exists with \c double parameters.
 *
 * @param fromPoints Points x as inhomogeneous 2-vectors
 * @param toPoints Points x' as inhomogeneous 2-vectors
 * @return calculated homography
 */
UBITRACK_EXPORT Math::Matrix< 3, 3, float > homographyDLT( const std::vector< Math::Vector< 2, float > >& fromPoints, 
	const std::vector< Math::Vector< 2, float > >& toPoints );

UBITRACK_EXPORT Math::Matrix< 3, 3, double > homographyDLT( const std::vector< Math::Vector< 2, double > >& fromPoints, 
	const std::vector< Math::Vector< 2, double > >& toPoints );
	

/**
 * @ingroup tracking_algorithms
 * Computes a homography for a square.
 *
 * The result is a homography H that maps vectors x = [(-0.5, +0.5), (-0.5, -0.5), (+0.5, -0.5), (+0.5, +0.5)] 
 * to the points x' given in the corners argument via x' = H x. The homography is computed
 * using the Harker & O'Leary method.
 * 
 * Note: also exists with \c double parameters.
 *
 * @param corners list of four 2D points x'
 * @return the homography
 */
UBITRACK_EXPORT Math::Matrix< 3, 3, float > squareHomography( const std::vector< Math::Vector< 2, float > >& corners );

UBITRACK_EXPORT Math::Matrix< 3, 3, double > squareHomography( const std::vector< Math::Vector< 2, double > >& corners );

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

#endif

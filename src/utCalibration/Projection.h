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
 * Functions dealing with projection matrices
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_CALIBRATION_PROJECTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_PROJECTION_H_INCLUDED__



#include <utCore.h>
#include <utMath/Matrix.h>
#include <utMath/Vector.h>
#include <vector>

namespace Ubitrack { namespace Calibration {

#ifdef HAVE_LAPACK

/**
 * @ingroup tracking_algorithms
 * Computes a 3x4 projection matrix using a linear DLT method.
 *
 * The result is a projection matrix P that maps 3d-points x to 2d-points x' via x'=Px.
 * See Hartley&Zisserman for details.
 *
 * Note: also exists with \c double parameters
 *
 * @param fromPoints Points x as inhomogeneous 3-vectors
 * @param toPoints Points x' as inhomogeneous 2-vectors
 * @return calculated projection matrix
 */
UBITRACK_EXPORT Math::Matrix< 3, 4, float > projectionDLT( const std::vector< Math::Vector< float, 3 > >& fromPoints, 
	const std::vector< Math::Vector< float, 2 > >& toPoints );

UBITRACK_EXPORT Math::Matrix< 3, 4, double > projectionDLT( const std::vector< Math::Vector< double, 3 > >& fromPoints, 
	const std::vector< Math::Vector< double, 2 > >& toPoints );


/**
 * @ingroup tracking_algorithms
 * Decomposes a 3x4 projection matrix into matrices K and R and a translation vector t.
 *
 * If you create a matrix/pose from r|t, the expression x' = [r|t] x will transform world 
 * coordinates x to camera coordinates x'. To get the camera position, invert this pose/matrix.
 *
 * Note: also exists with \c double parameters
 *
 * @param k resulting camera intrinsics matrix (upper triangular)
 * @param r resulting orthogonal rotation matrix
 * @param t resulting translation vector
 * @param p the input projection matrix
 */
UBITRACK_EXPORT void decomposeProjection( Math::Matrix< 3, 3, float >& k, 
	Math::Matrix< 3, 3, float >& r, Math::Vector< float, 3 >& t, const Math::Matrix< 3, 4, float >& p ); 

UBITRACK_EXPORT void decomposeProjection( Math::Matrix< 3, 3, double >& k, 
	Math::Matrix< 3, 3, double >& r, Math::Vector< double, 3 >& t, const Math::Matrix< 3, 4, double >& p );    

#endif // HAVE_LAPACK

/**
 * @ingroup tracking_algorithms
 * Computes a 4x4 projection matrix for OpenGL.
 *
 * You can use a 3x3 odr 4x3 projection matrix for the computation.
 *
 * @param l left border of the camera image
 * @param r right border of the camera image
 * @param b bottom border of the camera image
 * @param t top border of the camera image 
 * @param n near clipping plane
 * @param f far clipping plane
 * @param matrix the input projection matrix
 */
UBITRACK_EXPORT Math::Matrix< 4, 4, double > projectionMatrixToOpenGL( double l, double r, double b, double t, double n, double f, Math::Matrix< 3, 4, double > m );

UBITRACK_EXPORT Math::Matrix< 4, 4, float > projectionMatrixToOpenGL( float l, float r, float b, float t, float n, float f, Math::Matrix< 3, 4, float > m );

UBITRACK_EXPORT Math::Matrix< 4, 4, double > projectionMatrixToOpenGL( double l, double r, double b, double t, double n, double f, Math::Matrix< 3, 3, double > m );

UBITRACK_EXPORT Math::Matrix< 4, 4, float > projectionMatrixToOpenGL( float l, float r, float b, float t, float n, float f, Math::Matrix< 3, 3, float > m );

/**
 * @ingroup tracking_algorithms
 * Computes a 4x4 off-axis projection matrix for OpenGL.
 *
 * @param eye viewpoint position
 * @param ll lower left egde of screen
 * @param ul upper left edge of screen
 * @param lr lower right edge of screen
 * @param n near clipping plane
 * @param f far clipping plane
 * @param sw screen width
 * @param sh screen height
 */
UBITRACK_EXPORT Math::Matrix< 4, 4, double > offAxisProjectionMatrix( Math::Vector< double, 3 >& eye, Math::Vector< double, 3 >& ll, Math::Vector< double, 3 >& ul, Math::Vector< double, 3 >& lr, double n, double f, double sw, double sh );


/**
 * @ingroup tracking_algorithms
 * Corrects a 3x3 projection matrix according to the origin flag of an image.
 * Note: In Ubitrack, projection matrices assume a bottom-up orientation!
 * @param k the intrinsic matrix to correct
 * @param origin the IplImage origin flag
 * @param height of the image in pixels
 */
UBITRACK_EXPORT void correctOrigin( Math::Matrix< 3, 3, float >& k, int origin, int height );
 
UBITRACK_EXPORT void correctOrigin( Math::Matrix< 3, 3, double >& k, int origin, int height );

} } // namespace Ubitrack::Calibration

#endif

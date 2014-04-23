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
 * Functions for 3D Point Reconstruction
 *
 * @author Daniel Muhra <muhra@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_CALIBRATION_3DPOINTRECONSTRUCTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_3DPOINTRECONSTRUCTION_H_INCLUDED__


#include <utCore.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>

namespace Ubitrack { namespace Algorithm {

/**
 * @ingroup tracking_algorithms
 * Computes the distance between a point and the epipole of the other point in the same picture
 *
 * The result is a relative distance of the two points
 *
 * Note: also exists with \c double parameters. You can either use 2D coordinates or homogeneous coordinates ( 3D )
 *
 * @param fromPoints Points x as inhomogeneous 2-vectors
 * @param toPoints Points x' as inhomogeneous 2-vectors
 * @param stepSize optional parameter, if you want to use e.g. only every second value for computation, then set stepSize to 2
 * @return calculated fundamental matrix
 */
UBITRACK_EXPORT float pointToPointDist( const Math::Vector< float, 2 > & from, const Math::Vector< float, 2 > & to, const Math::Matrix< float, 3, 3 > & fM  );

UBITRACK_EXPORT double pointToPointDist( const Math::Vector< double, 2 > & from, const Math::Vector< double, 2 > & to, const Math::Matrix< double, 3, 3 > & fM  );

UBITRACK_EXPORT float pointToPointDist( const Math::Vector< float, 3 > & from, const Math::Vector< float, 3 > & to, const Math::Matrix< float, 3, 3 > & fM  );

UBITRACK_EXPORT double pointToPointDist( const Math::Vector< double, 3 > & from, const Math::Vector< double, 3 > & to, const Math::Matrix< double, 3, 3 > & fM  );

#ifdef HAVE_LAPACK
/**
 * @ingroup tracking_algorithms
 * Estimates the 3D position of a point seen by two cameras
 *
 * The result is a 3D vector( position ) of the point
 *
 * Note: also exists with \c double parameters. You can either use 2D coordinates or homogenous coordinates ( 3D )
 *
 * @param P1 the projection matrix of the first camera
 * @param P2 the projection matrix of the second camera
 * @param p1 the position of the point in the image of the first camera
 * @param p2 the position of the point in the image of the second camera
 * @return estimated 3D position
 */

UBITRACK_EXPORT Math::Vector< float, 3 > get3DPosition( const Math::Matrix< float, 3, 4 > & P1, const Math::Matrix< float, 3, 4 > & P2, const Math::Vector< float, 2 > & p1, const Math::Vector< float, 2 > & p2 );

UBITRACK_EXPORT Math::Vector< double, 3 > get3DPosition( const Math::Matrix< double, 3, 4 > & P1, const Math::Matrix< double, 3, 4 > & P2, const Math::Vector< double, 2 > & p1, const Math::Vector< double, 2 > & p2 );

/**
 * @ingroup tracking_algorithms
 * Reconstructs 3D points from two sets of 2D points
 *
 * The result is a vector of 3D points
 *
 * Note: also exists with \c double parameters. You can either use 2D coordinates or homogenous coordinates ( 3D )
 *
 * @param p1 a vector of 2D points from the first camera
 * @param p2 a vector of 2D points from the second camera
 * @param P1 the projection matrix of the first camera
 * @param P2 the projection matrix of the second camera
 * @param F the fundamental matrix
 * @return a vector of 3D points
 */
UBITRACK_EXPORT std::vector< Math::Vector< float, 3 > > reconstruct3DPoints( const std::vector< Math::Vector< float, 2 > > & p1, const std::vector< Math::Vector< float, 2 > > & p2,
																			const Math::Matrix< float, 3, 4 > & P1, const Math::Matrix< float, 3, 4 > & P2, const Math::Matrix< float, 3, 3 > & fM );

UBITRACK_EXPORT std::vector< Math::Vector< double, 3 > > reconstruct3DPoints( const std::vector< Math::Vector< double, 2 > > & p1, const std::vector< Math::Vector< double, 2 > > & p2,
																			const Math::Matrix< double, 3, 4 > & P1, const Math::Matrix< double, 3, 4 > & P2, const Math::Matrix< double, 3, 3 > & fM );


/**
 * @ingroup tracking_algorithms
 * Reconstructs a 3D point from n camera matrices ( n > 1 ) and n observations
 *
 * The result is the 3D position
 *
 * Note: also exists with \c double parameters.
 *
 * @param P a vector of camera matrices
 * @param points a vector of 2D points from the cameras
 * @param falg to indicate minimization
 * flag == 0: algebraic(fast), flag == 1: non-linear (accurate)
 * @return a 3D points
 */

UBITRACK_EXPORT Math::Vector< float, 3 > get3DPosition( const std::vector< Math::Matrix< float, 3, 4 > > &P, const std::vector< Math::Vector< float, 2 > > &points, std::size_t flag );

UBITRACK_EXPORT Math::Vector< double, 3 > get3DPosition( const std::vector< Math::Matrix< double, 3, 4 > > &P, const std::vector< Math::Vector< double, 2 > >& points, std::size_t flag );

UBITRACK_EXPORT Math::Vector< double, 3 > get3DPositionWithResidual( const std::vector< Math::Matrix< double, 3, 4 > > &P, const std::vector< Math::Vector< double, 2 > >& points, std::size_t flag = 0, double* residual = 0 );
#endif

} } // namespace Ubitrack::Algorithm

#endif

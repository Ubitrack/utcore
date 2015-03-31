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
 * cpp-file for lens un/-distortion.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 
 

#include "Correction.h"
#include "Distortion.h"

#ifdef HAVE_LAPACK
#include "Undistortion.h"
#endif

namespace Ubitrack { namespace Algorithm { namespace CameraLens {


void distort( const Ubitrack::Math::CameraIntrinsics< float >& mat, const Math::Vector2f& undistorted, Math::Vector2f& distorted )
{
	distort_impl( mat, undistorted, distorted );
}

void distort( const Ubitrack::Math::CameraIntrinsics< double >& mat, const Math::Vector2d& undistorted, Math::Vector2d& distorted )
{
	distort_impl( mat, undistorted, distorted );	
}

void distort( const Ubitrack::Math::CameraIntrinsics< float >& mat, const std::vector< Math::Vector2f >& undistorted, std::vector< Math::Vector2f >& distorted )
{
	distort_impl( mat, undistorted, distorted );	
}

void distort( const Ubitrack::Math::CameraIntrinsics< double >& mat, const std::vector< Math::Vector2d >& undistorted, std::vector< Math::Vector2d >& distorted )
{
	distort_impl( mat, undistorted, distorted );	
}
 

#ifdef HAVE_LAPACK


void undistort( const Ubitrack::Math::CameraIntrinsics< float >& mat, const Math::Vector2f& distorted, Math::Vector2f& undistorted )
{
	undistort_impl( mat, distorted, undistorted );
}

void undistort( const Ubitrack::Math::CameraIntrinsics< double >& mat, const Math::Vector2d& distorted, Math::Vector2d& undistorted )
{
	undistort_impl( mat, distorted, undistorted );
}

void undistort( const Ubitrack::Math::CameraIntrinsics< float >& mat, const std::vector< Math::Vector2f >& distorted, std::vector< Math::Vector2f >& undistorted )
{
	undistort_impl( mat, distorted, undistorted );	
}

void undistort( const Ubitrack::Math::CameraIntrinsics< double >& mat, const std::vector< Math::Vector2d >& distorted, std::vector< Math::Vector2d >& undistorted )
{
	undistort_impl( mat, distorted, undistorted );
}
	
#endif
	
}}} // namespace Ubitrack::Algorithm::CameraLens
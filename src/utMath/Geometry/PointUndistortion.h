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
 *
 * Several functions that can be applied to points or container of points.
 *
 * The header includes functions for common operations on points.
 * The functions can easily be combined to STL-containers like 
 * \c std::vectors, \c std::list, etc, that contain points.
 *
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __UBITRACK_MATH_GEOMETRY_POINT_UNDISTORTION_H__
#define __UBITRACK_MATH_GEOMETRY_POINT_UNDISTORTION_H__

// Ubitrack
#include <utMath/Vector.h>
#include <utMath/CameraIntrinsics.h>

namespace Ubitrack { namespace Math { namespace Geometry {

template < typename T >
struct PointUndistortion
{

	const Math::CameraIntrinsics< T > &intrinsics;

public:
	PointUndistortion( const Math::CameraIntrinsics< T >& _intrinsics )
		: intrinsics( _intrinsics )
		{
		}
		
	// internal template to catch wrong vector types and print an error message
	template< typename notSupportedVectorType >
	notSupportedVectorType operator() ( const notSupportedVectorType& vec_in ) const
	{
		UBITRACK_STATIC_ASSERT( false, USE_ONLY_WITH_VECTOR_TYPE_OF_2_DIMENSIONS );
		return vec_in;
	}
	
	
	/// @internal Specialization of \c bracket-operator for \b 2D \b point undistortion
	template< typename VType >
	void operator() ( const Math::Vector< VType, 2 > &vecIn, Math::Vector< VType, 2 > &vecOut ) const
	{
		// 1.) take images to sensor coordinate system
		// [x y]^T = [K]^-1 * [x' y' 1]^T
		const VType x( ( ( vecIn[ 0 ] + intrinsics.matrix( 0, 2 ) ) / intrinsics.matrix( 0, 0 ) ) );
		const VType y( ( ( vecIn[ 1 ] + intrinsics.matrix( 1, 2 ) ) / intrinsics.matrix( 1, 1 ) ) );
		
		// 2.) generate the polynomial factors for distortion, depending on the radius 
		// r^2=(x^2+y^2), r^4, r^6
		const VType radiusPow2 = ( x*x + y*y );
		const VType radiusPow4 = radiusPow2*radiusPow2;
		const VType radiusPow6 = radiusPow4*radiusPow2;
		
		// 3.) calculate the scaling factor of the radial distortion
		//	f=(1+r^2k_1+1+r^4k_2+1+r^6k_3 )/ (1+r^2k_4+1+r^4k_5+1+r^6k_6 )
		const VType nom = (1+ intrinsics.radial_params[ 0 ] * radiusPow2 + intrinsics.radial_params[ 1 ]  * radiusPow4 + intrinsics.radial_params[ 2 ]  * radiusPow6 );
		const VType denom = (1+ intrinsics.radial_params[ 3 ] * radiusPow2 + intrinsics.radial_params[ 4 ]  * radiusPow4 + intrinsics.radial_params[ 5 ]  * radiusPow6 );
		const VType scale = nom/denom;
		
		// 4.) scale points by distortion factor and move on tangential line as well
		// x'' = x*f + 2*t_1*xy + t_2*(r^2+2x^2)
		// y'' = y*f + 2*t_2*xy + t_1*(r^2+2y^2)
		const VType x_undis = x * scale + 2 * intrinsics.tangential_params[ 0 ] * x * y + intrinsics.tangential_params[ 1 ] *( radiusPow2+2*x*x );
		const VType y_undis = y * scale + 2 * intrinsics.tangential_params[ 1 ] * x * y + intrinsics.tangential_params[ 0 ] *( radiusPow2+2*y*y );
		
		// 5.) project point into image plane again
		// [u,v]^T = [K] * [x'' y'' 1]^T
		vecOut[ 0 ] = ( intrinsics.matrix( 0, 0 ) * x_undis - intrinsics.matrix( 0, 2 ) );
		vecOut[ 1 ] = ( intrinsics.matrix( 1, 1 ) * y_undis - intrinsics.matrix( 1, 2 ) );
		
	}
	
	template< typename VType >
	Math::Vector< VType, 2 > operator() ( const Math::Vector< VType, 2 > &vecIn ) const
	{
		Math::Vector< T, 2 > vecOut;
		this->operator()( vecIn, vecOut );
		return vecOut;
	}
	
};

/// function wrapper that encapsulates the call to the functor @todo more documentation here
template< typename T, typename ForwardIterator1, typename ForwardIterator2 >
inline void undistort_points( const Math::CameraIntrinsics< T > &intrinsics, const ForwardIterator1 iBegin, const ForwardIterator1 iEnd, ForwardIterator2 iOut )
{
	std::transform( iBegin, iEnd, iOut, PointUndistortion< T >( intrinsics ) );
}

/// function wrapper that encapsulates the call to the functor @todo more documentation here
template< typename T, typename VType >
inline void undistort_point( const Math::CameraIntrinsics< T >& camIntrin, const Math::Vector< VType, 2 > &vecIn, Math::Vector< VType, 2 > &vecOut )
{
	const PointUndistortion< T >tmp( camIntrin );//( vecIn, vecOut );
	tmp( vecIn, vecOut );
}

}}} // namespace Ubitrack::Math::Geometry

#endif // __UBITRACK_MATH_GEOMETRY_POINT_UNDISTORTION_H__
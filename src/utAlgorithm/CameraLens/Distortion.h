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
 * Implementation of radial (and tangential) 2-vector(s) distortion
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_CAMERALENS_DISTORTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_CAMERALENS_DISTORTION_H_INCLUDED__


namespace Ubitrack { namespace Algorithm { namespace CameraLens {
	
namespace internal {
	
	/// @internal helper function to unproject a single point from image(pixel) to sensor coordinates)
	template< typename T, std::size_t N >
	inline void unproject_impl( const Math::CameraIntrinsics< T >& camIntrin, const Math::Vector< T, N >& imagePt, Math::Vector< T, N >& result )
	{
		using namespace Ubitrack;
		const Math::Matrix< T, 3, 3 >& camMat = camIntrin.matrix;
		const Math::Vector< T, 2 >& pt = imagePt;
		// unproject point to normalized camera coordinates (assuming K( 2, 2 ) == +/-1)
		T y = ( pt( 1 ) - camMat( 1, 2 ) * camMat( 2, 2 ) ) / camMat( 1, 1 );
		T x = ( pt( 0 ) - camMat( 0, 1 ) * y - camMat( 0, 2 ) * camMat( 2, 2 ) ) / camMat( 0, 0 );
		result( 0 ) = x;
		result( 1 ) = y;
	}

	/// @internal helper function to project a single point from sensor to image(pixel) coordinates)
	template< typename T, std::size_t N >
	inline void project_impl( const Math::CameraIntrinsics< T >& camIntrin, const Math::Vector< T, N >& sensorPt, Math::Vector< T, N >& result )
	{
		const Math::Matrix< T, 3, 3 >& camMat = camIntrin.matrix;
		const Math::Vector< T, 2 >& pt = sensorPt;
		const T x = pt( 0 ) * camMat( 0, 0 ) + pt( 1 ) * camMat( 0, 1 ) + camMat( 0, 2 ) * camMat( 2, 2 );
		const T y = pt( 1 ) * camMat( 1, 1 ) + camMat( 1, 2 ) * camMat( 2, 2 );
		result( 0 ) = x;
		result( 1 ) = y;
	}

	/// @internal implements the distortion function for general types
	template< class VT1, class VT2, class VT3, class VT4 >
	inline void distort_impl( const VT1& radVec, const VT2& tanVec, const VT3& pointIn, VT4& result )
	{
		typedef typename VT1::value_type VType;
		const VType x = pointIn( 0 );
		const VType y = pointIn( 1 );
		const VType xx = x*x;
		const VType xy2 = x*y*2;
		const VType yy = y*y;
		const VType r2 = xx+yy;
		const VType r4 = r2*r2;
		const VType r6 = r4*r2;
		const VType upper = ( 1 + radVec( 0 ) * r2 + radVec( 1 ) * r4 + radVec( 2 ) * r6 );
		const VType lower = ( 1 + radVec( 3 ) * r2 + radVec( 4 ) * r4 + radVec( 5 ) * r6 );
		const VType ratio = upper/lower;
		
		result( 0 ) = (x * ratio) +  tanVec ( 0 ) * xy2 + tanVec( 1 ) * ( r2 + xx * 2 );
		result( 1 ) = (y * ratio) +  tanVec ( 1 ) * xy2 + tanVec( 0 ) * ( r2 + yy * 2 );
	}


	template< typename T, std::size_t N >
	inline void distort_impl( const Math::CameraIntrinsics< T >& camIntrin,  const Math::Vector< T, N >& undistorted, Math::Vector< T, N >& distorted )
	{
		Math::Vector< T, N > result;
		unproject_impl( camIntrin, undistorted, result );
		
		distort_impl( camIntrin.radial_params, camIntrin.tangential_params, result, distorted );
		project_impl( camIntrin, distorted, distorted );
	}

}	// namespsce ::internal 

template< typename T, std::size_t N >
inline void distort( const Math::Vector< T, 6 >& radVec, const Math::Vector< T, 2 >& tanVec, const Math::Vector< T, N >& pointIn,  Math::Vector< T, N >& result )
{
	internal::distort_impl( radVec, tanVec, pointIn, result );
}

template< typename RadVec, typename TanVec, typename UndisType, typename DisType >
inline void distort( const RadVec& radVec, const TanVec& tanVec, const UndisType& pointIn, DisType& result )
{
	internal::distort_impl( radVec, tanVec, pointIn, result );
}

template< typename T, std::size_t N >
inline void distort( const Math::Vector< T, 6 >& radVec, const Math::Vector< T, 2 >& tanVec, const std::vector< Math::Vector< T, N > >& pointsIn, std::vector< Math::Vector< T, N > >& result )
{
	typename std::vector< Math::Vector< T, N > >::const_iterator itBegin = pointsIn.begin();
	const typename std::vector< Math::Vector< T, N > >::const_iterator itEnd = pointsIn.end();
	/// @todo check here if same vector is used as output or another -> could be std::back_inserter than
	typename std::vector< Math::Vector< T, N > >::iterator itOut = result.begin();
	for( ; itBegin != itEnd; ++itBegin, ++itOut )
		internal::distort_impl( radVec, tanVec, *itBegin, *itOut );
}

template< typename T, std::size_t N >
inline void distort_impl( const Math::CameraIntrinsics< T >& camIntrin,  const Math::Vector< T, N >& undistorted, Math::Vector< T, N >& distorted )
{
	internal::distort_impl( camIntrin, undistorted, distorted );
}


template< typename T, std::size_t N >
inline void distort_impl( const Math::CameraIntrinsics< T > camIntrin, const std::vector< Math::Vector< T, N > >& pointsIn, std::vector< Math::Vector< T, N > >& result )
{
	typename std::vector< Math::Vector< T, N > >::const_iterator itBegin = pointsIn.begin();
	const typename std::vector< Math::Vector< T, N > >::const_iterator itEnd = pointsIn.end();
	/// @todo check here if same vector is used as output or another -> could be std::back_inserter than
	typename std::vector< Math::Vector< T, N > >::iterator itOut = result.begin();
	for( ; itBegin != itEnd; ++itBegin, ++itOut )
		internal::distort_impl( camIntrin, *itBegin, *itOut );
}


}}} // namespace Ubitrack::Algorithm::CameraLens

#endif //__UBITRACK_CALIBRATION_FUNCTION_CAMERALENS_DISTORTION_H_INCLUDED__

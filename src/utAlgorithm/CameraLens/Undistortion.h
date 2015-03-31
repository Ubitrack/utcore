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
 * Implementation of radial (and tangential) 2-vector(s) undistortion
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_CAMERALENS_UNDISTORTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_CAMERALENS_UNDISTORTION_H_INCLUDED__

#include "Distortion.h"
#include <utMath/CameraIntrinsics.h>

#ifdef HAVE_LAPACK
#include <utMath/Optimization/LevenbergMarquardt.h>
#endif

namespace Ubitrack { namespace Algorithm { namespace CameraLens {

namespace internal {
/**
 * Radially and tangentially undistorts a 2-vector (x,y). The distortion is described by two vectors [k1, k2, k3, k4, k5, k6 ] and [ p1, p2].
 *    x' = x( ( 1 + k1 * r^2 + k2 * r^4 + k3 * r^6 ) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6) ) + ( 2 * p1 * x * y + p2 * ( r^2 + 2 * x^2 ) )
 *    y' = x( ( 1 + k1 * r^2 + k2 * r^4 + k3 * r^6 ) / (1 + k4 * r^2 + k5 * r^4 + k6 * r^6) ) + ( 2 * p2 * x * y + p1 * ( r^2 + 2 * y^2 ) )
 * where
 *    r^2 = x^2 + y^2
 *
 *
 * The jacobian is computed wrt. the point p.
 */
template< typename T >
class PointUndistortion
{
protected:
	const Math::Vector< T, 6 >& m_k;
	const Math::Vector< T, 2 >& m_p;
	
public:
	/**
	 * Constructor.
	 * @param d reference to the distortion coefficients
	 */
	 
	PointUndistortion( const Math::CameraIntrinsics< T >& cam )
		: m_k( cam.radial_params )
		, m_p( cam.tangential.radial_params )
	{}
	
	PointUndistortion( const Math::Vector< T, 6 >& radVec, const Math::Vector< T, 2 >& tanVec )
		: m_k( radVec )
		, m_p( tanVec )
	{}

	/**
	 * return the size of the result vector
	 */
	std::size_t size() const
	{
		return 2;
	}

	/**
	 * @param result the transformed point (2-vector)
	 * @param input the 2-vector describing the undistorted point
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		distort( m_k, m_p, input, result );
	}
	
	/**
	 * @param result the transformed point (2-vector)
	 * @param input the 2-vector describing the undistorted point
	 * @param J jacobian (2x2-matrix)
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input the 2-vector describing the undistorted point
	 * @param J jacobian (2x2-matrix)
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		// the following matlab code (symbolic toolbox) was used to generate the derivative
		// syms x y k1 k2 k2 k3 k4 k5 k6 p1 p2
		// r2 = x*x + y*y
		// % xx = x + x * ( k1 * r2 + k2 * r2*r2 ) + ( 2 * p1 * x * y + p2 * ( r2 + 2 * x*x ) )
		// % yy = y + y * ( k1 * r2 + k2 * r2*r2 ) + ( 2 * p2 * x * y + p1 * ( r2 + 2 * y*y ) )
		// % xx = x + x * ( k1 * r2 + k2 * r2*r2 + k3*r2*r2*r2) + ( 2 * p1 * x * y + p2 * ( r2 + 2 * x*x ) ) % for 3 params
		// % yy = y + y * ( k1 * r2 + k2 * r2*r2 + k3*r2*r2*r2) + ( 2 * p2 * x * y + p1 * ( r2 + 2 * y*y ) ) % for 3 params
		// xx = x * (( 1 + k1 * r2 + k2 * r2*r2 + k3*r2*r2*r2) / ( 1+ k4 * r2 + k5 * r2*r2 + k6*r2*r2*r2) ) + ( 2 * p1 * x * y + p2 * ( r2 + 2 * x*x ) ) % for 6 params
		// yy = y * (( 1 + k1 * r2 + k2 * r2*r2 + k3*r2*r2*r2) / ( 1+ k4 * r2 + k5 * r2*r2 + k6*r2*r2*r2) ) + ( 2 * p2 * x * y + p1 * ( r2 + 2 * y*y ) ) % for 6 params
		// ccode( jacobian( [xx; yy], [x,y] ) )
		// remark: the 3rd radial distortion parameter is integrated here, when it is null the term will have no effect
		// the same for the other 3 (later) parameters
		
		typedef typename MT::value_type VType;
		const VType x = input( 0 );
		const VType y = input( 1 );
		const VType k1 = m_k( 0 );
		const VType k2 = m_k( 1 );
		const VType k3 = m_k( 2 );
		const VType k4 = m_k( 3 );
		const VType k5 = m_k( 4 );
		const VType k6 = m_k( 5 );
		const VType p1 = m_p( 0 );
		const VType p2 = m_p( 1 );
		
		const VType t3 = x*x;
		const VType t4 = y*y;
		const VType t2 = t3+t4;
		const VType t5 = t2*t2;
		const VType t6 = k5*t5;
		const VType t7 = k6*t2*t5;
		const VType t8 = k4*t2;
		const VType t9 = t6+t7+t8+1;
		const VType t10 = 1/t9;
		const VType t11 = k2*t5;
		const VType t12 = k3*t2*t5;
		const VType t13 = k1*t2;
		const VType t14 = t11+t12+t13+1;
		const VType t15 = 1/(t9*t9);
		const VType t16 = p1*x*2;
		const VType t17 = p2*y*2;
		const VType t18 = k1*x*2;
		const VType t19 = k2*t2*x*4;
		const VType t20 = k3*t5*x*6;
		const VType t21 = t18+t19+t20;
		const VType t22 = k4*x*2;
		const VType t23 = k5*t2*x*4;
		const VType t24 = k6*t5*x*6;
		const VType t25 = t22+t23+t24;
		const VType t26 = t10*t14;
		const VType t27 = k1*y*2;
		const VType t28 = k2*t2*y*4;
		const VType t29 = k3*t5*y*6;
		const VType t30 = t27+t28+t29;
		const VType t31 = k4*y*2;
		const VType t32 = k5*t2*y*4;
		const VType t33 = k6*t5*y*6;
		const VType t34 = t31+t32+t33;
		J( 0, 0 ) = t26+p2*x*6+p1*y*2+t10*t21*x-t14*t15*t25*x;
		J( 0, 1 ) = t16+t17+t10*t30*x-t14*t15*t34*x;
		J( 1, 0 ) = t16+t17+t10*t21*y-t14*t15*t25*y;
		J( 1, 1 ) = t26+p2*x*2+p1*y*6+t10*t30*y-t14*t15*t34*y;
		
		/// for 2 radial distortion coefficients
		// J( 0, 0 ) = k2*pow(x*x+y*y,2.0)+p2*x*6.0+p1*y*2.0+x*(k1*x*2.0+k2*x*(x*x+y*y)*4.0)+k1*(x*x+y*y)+1.0;
		// J( 0, 1 ) = p1*x*2.0+p2*y*2.0+x*(k1*y*2.0+k2*y*(x*x+y*y)*4.0);
		// J( 1, 0 ) = p1*x*2.0+p2*y*2.0+y*(k1*x*2.0+k2*x*(x*x+y*y)*4.0);
		// J( 1, 1 ) = k2*pow(x*x+y*y,2.0)+p2*x*2.0+p1*y*6.0+y*(k1*y*2.0+k2*y*(x*x+y*y)*4.0)+k1*(x*x+y*y)+1.0;
		
		/// for 3 radial distortion coefficients
		// J( 0, 0 ) = k2*pow(x*x+y*y,2.0)+k3*pow(x*x+y*y,3.0)+p2*x*6.0+p1*y*2.0+x*(k1*x*2.0+k2*x*(x*x+y*y)*4.0+k3*x*pow(x*x+y*y,2.0)*6.0)+k1*(x*x+y*y)+1.0;
		// J( 0, 1 ) = p1*x*2.0+p2*y*2.0+x*(k1*y*2.0+k2*y*(x*x+y*y)*4.0+k3*y*pow(x*x+y*y,2.0)*6.0);
		// J( 1, 0 ) = p1*x*2.0+p2*y*2.0+y*(k1*x*2.0+k2*x*(x*x+y*y)*4.0+k3*x*pow(x*x+y*y,2.0)*6.0);
		// J( 1, 1 ) = k2*pow(x*x+y*y,2.0)+k3*pow(x*x+y*y,3.0)+p2*x*2.0+p1*y*6.0+y*(k1*y*2.0+k2*y*(x*x+y*y)*4.0+k3*y*pow(x*x+y*y,2.0)*6.0)+k1*(x*x+y*y)+1.0;

		/// for six radial distortion coefficients
		// J( 0, 0 ) = (k2*pow(x*x+y*y,2)+k3*pow(x*x+y*y,3)+k1*(x*x+y*y)+1)/(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1)+p2*x*6+p1*y*2+(x*(k1*x*2+k2*x*(x*x+y*y)*4+k3*x*pow(x*x+y*y,2)*6))/(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1)-x*(k4*x*2+k5*x*(x*x+y*y)*4+k6*x*pow(x*x+y*y,2)*6)*(k2*pow(x*x+y*y,2)+k3*pow(x*x+y*y,3)+k1*(x*x+y*y)+1)*1/pow(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1,2);
		// J( 0, 1 ) = p1*x*2+p2*y*2+(x*(k1*y*2+k2*y*(x*x+y*y)*4+k3*y*pow(x*x+y*y,2)*6))/(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1)-x*(k4*y*2+k5*y*(x*x+y*y)*4+k6*y*pow(x*x+y*y,2)*6)*(k2*pow(x*x+y*y,2)+k3*pow(x*x+y*y,3)+k1*(x*x+y*y)+1)*1/pow(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1,2);
		// J( 1, 0 ) = p1*x*2+p2*y*2+(y*(k1*x*2+k2*x*(x*x+y*y)*4+k3*x*pow(x*x+y*y,2)*6))/(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1)-y*(k4*x*2+k5*x*(x*x+y*y)*4+k6*x*pow(x*x+y*y,2)*6)*(k2*pow(x*x+y*y,2)+k3*pow(x*x+y*y,3)+k1*(x*x+y*y)+1)*1/pow(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1,2);
		// J( 1, 1 ) = (k2*pow(x*x+y*y,2)+k3*pow(x*x+y*y,3)+k1*(x*x+y*y)+1)/(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1)+p2*x*2+p1*y*6+(y*(k1*y*2+k2*y*(x*x+y*y)*4+k3*y*pow(x*x+y*y,2)*6))/(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1)-y*(k4*y*2+k5*y*(x*x+y*y)*4+k6*y*pow(x*x+y*y,2)*6)*(k2*pow(x*x+y*y,2)+k3*pow(x*x+y*y,3)+k1*(x*x+y*y)+1)*1/pow(k5*pow(x*x+y*y,2)+k6*pow(x*x+y*y,3)+k4*(x*x+y*y)+1,2);
	}
};
}	// namespace ::internal

#ifdef HAVE_LAPACK

template< typename T, std::size_t N >
void undistort_impl( const Math::Vector< T, 6 >& radVector, const Math::Vector< T, 2 >& tanVector, const Math::Vector< T, N >& distorted, Math::Vector< T, N >& undistorted )
{
	// non-linear minimization
	Math::Vector< T, N > tmp = ( distorted );
	internal::PointUndistortion< T > distFunc( radVector, tanVector );
	Math::Optimization::levenbergMarquardt( distFunc, tmp, distorted, Math::Optimization::OptTerminate( 5, 1e-5 ), Math::Optimization::OptNoNormalize() );
	undistorted = tmp;
}

template< typename T >
void undistort_impl( const Math::CameraIntrinsics< T >& mat, const Math::Vector< T, 2 >& distorted, Math::Vector< T, 2 >& undistorted )
{
	using namespace Ubitrack;
	
	Math::Vector< T, 2 > camPoint;
	internal::unproject_impl( mat, distorted, camPoint );
	// non-linear minimization
	undistort_impl( mat.radial_params, mat.tangential_params, camPoint, undistorted );
	
	// back to image coordinates
	internal::project_impl( mat, undistorted, undistorted );
}

template< typename T, std::size_t N >
inline void undistort_impl( const Math::CameraIntrinsics< T > camIntrin, const std::vector< Math::Vector< T, N > >& pointsIn, std::vector< Math::Vector< T, N > >& result )
{
	typename std::vector< Math::Vector< T, N > >::const_iterator itBegin = pointsIn.begin();
	const typename std::vector< Math::Vector< T, N > >::const_iterator itEnd = pointsIn.end();
	/// @todo check here if same vector is used as output or another -> could be std::back_inserter than
	typename std::vector< Math::Vector< T, N > >::iterator itOut = result.begin();
	for( ; itBegin != itEnd; ++itBegin, ++itOut )
		undistort_impl( camIntrin, *itBegin, *itOut );
}

#endif

}}} // namespace Ubitrack::Algorithm::CameraLens

#endif //__UBITRACK_CALIBRATION_FUNCTION_CAMERALENS_UNDISTORTION_H_INCLUDED__

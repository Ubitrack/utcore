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
 * Functions for ransac absolute orientation estimation .
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_OPTIMIZATION_H_INCLUDED__
#define __UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_OPTIMIZATION_H_INCLUDED__


#include <utCore.h>
#include <utMath/Pose.h>
#include <utMath/Blas1.h>
#include <utMath/Optimization/Optimization.h>

#ifdef HAVE_LAPACK
#include <utMath/Optimization/LevenbergMarquardt.h>
#endif	// HAVE_LAPACK

// Boost
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation3D3D {


/**
 * @internal 
 * A minimization function for the non-linear optimization
 * can be applied with Levenberg-Marquardt
 */
template< typename InputIterator >
class PointCorrespodencesSinglePose
{

	typedef typename std::iterator_traits< InputIterator >::value_type::value_type VType;
	
protected:
	const InputIterator m_iBegin;
	const InputIterator m_iEnd;
	
public:
	/** 
	 * constructor.
	 * @param iBegin iterator to the beginning of a container with 3d-points(must stay constant during lifetime of the object)
	 * @param iEnd iterator to the end of a container with 3d-point(must stay constant during lifetime of the object)
	 */
	PointCorrespodencesSinglePose( const InputIterator iBegin, const InputIterator iEnd )
		: m_iBegin( iBegin )
		, m_iEnd( iEnd )
	{}

	/**
	 * return the size of the result vector
	 */
	std::size_t size() const
	{ return ( 3 * std::distance( m_iBegin, m_iEnd ) ); }

	/**
	 * @param result 6-vector to store the result in
	 * @param input containing the pose parameters from one point cloud to the other one
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		// const Math::Pose pose( Math::Pose::fromVector( input ) );
		const Math::Pose pose = Math::Pose( 
		Math::Quaternion( input[ 3 ], input[ 4 ], input[ 5 ], input[ 6 ] ).normalize()
		, Math::Vector< VType, 3 >( input [ 0 ], input [ 1 ], input [ 2 ] )
		);
		
		
		std::size_t i ( 0 );
		for ( InputIterator it ( m_iBegin ); it != m_iEnd; ++i, ++it )
		{
			const Math::Vector< VType, 3 > pt( (*it)[ 0 ], (*it)[ 1 ], (*it)[ 2 ] );
			const Math::Vector< VType, 3 > vec = pose * pt;
			result[ i*3+0 ] = vec[ 0 ];
			result[ i*3+1 ] = vec[ 1 ];
			result[ i*3+2 ] = vec[ 2 ];
		}
		// std::cout << "Result : \n" << result << std::endl;
	}
	
	/**
	 * @param result vector to store the result in
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the Jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		// TODO: implement as one function (more efficient)
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the Jacobian (evaluated for input) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		// expected order: tx, ty, tz, qx, qy, qz, qw.
		
		const VType tx = input[ 0 ];
		const VType ty = input[ 1 ];
		const VType tz = input[ 2 ];
		const VType n = std::sqrt( input[ 3 ]*input[ 3 ]+ input[ 4 ]*input[ 4 ] + input[ 5 ]*input[ 5 ] +input[ 6 ]*input[ 6 ] );
		const VType qx = input[ 3 ] / n;
		const VType qy = input[ 4 ] / n;
		const VType qz = input[ 5 ] / n;
		const VType qw = input[ 6 ] / n;
		
		
		std::size_t i( 0 );
		for ( InputIterator it ( m_iBegin ); it != m_iEnd; ++i, ++it )
		{
			const VType x = (*it)[ 0 ];
			const VType y = (*it)[ 1 ];
			const VType z = (*it)[ 2 ];
			
			const VType t2 = qx*y*2;
			const VType t3 = qw*z*2;
			const VType t4 = qx*x*2;
			const VType t5 = qy*y*2;
			const VType t6 = qz*z*2;
			const VType t7 = t4+t5+t6;
			const VType t8 = qw*x*2;
			const VType t9 = qy*z*2;
			const VType t15 = qz*y*2;
			const VType t10 = t8+t9-t15;
			const VType t11 = qx*z*2;
			const VType t12 = qz*x*2;
			const VType t13 = qw*y*2;
			const VType t14 = -t11+t12+t13;
			const VType t16 = qy*x*2;
			J( i*3+0, 0 ) = 1;
			J( i*3+0, 3 ) = t7;
			J( i*3+0, 4 ) = t2+t3-qy*x*2;
			J( i*3+0, 5 ) = t11-qz*x*2-qw*y*2;
			J( i*3+0, 6 ) = t10;
			J( i*3+1, 1 ) = 1;
			J( i*3+1, 3 ) = -t2-t3+t16;
			J( i*3+1, 4 ) = t7;
			J( i*3+1, 5 ) = t10;
			J( i*3+1, 6 ) = t14;
			J( i*3+2, 2 ) = 1;
			J( i*3+2, 3 ) = t14;
			J( i*3+2, 4 ) = -t8+t9+t15;
			J( i*3+2, 5 ) = t4+t5-t6;
			J( i*3+2, 6 ) = t2+t3-t16;
			
			J( i*3+0, 1 ) = J( i*3+0, 2 ) = 0;
			J( i*3+1, 0 ) = J( i*3+1, 2 ) = 0;
			J( i*3+2, 0 ) = J( i*3+2, 1 ) = 0;
			
			// const VType t2 = qy*qy;
			// const VType t3 = qz*qz;
			// const VType t4 = qw*qw;
			// const VType t5 = qx*qx;
			// const VType t6 = t2+t3+t4+t5;
			// const VType t7 = 1.0/(t6*t6);
			// const VType t8 = qw*t4*z;
			// const VType t9 = qw*t5*z;
			// const VType t10 = qw*t3*z;
			// const VType t11 = qz*t3*z;
			// const VType t12 = qz*t4*z;
			// const VType t13 = qz*t2*z;
			// const VType t14 = qy*t2*z;
			// const VType t15 = qy*t5*z;
			// const VType t16 = qy*t3*z;
			// const VType t17 = qx*t3*z;
			// const VType t18 = qw*qy*qz*z*2.0;
			// const VType t19 = qw*t4*y;
			// const VType t20 = qz*t4*x;
			// const VType t21 = qw*t5*y;
			// const VType t22 = qw*t2*y;
			// const VType t23 = qw*qx*qy*x*2.0;
			// const VType t24 = qx*qy*qz*y*2.0;
			// const VType t25 = qw*t4*x;
			// const VType t26 = qw*t5*x;
			// const VType t27 = qw*t2*x;
			// const VType t28 = qz*t4*y;
			// const VType t29 = qy*t4*z;
			// const VType t30 = qx*t5*x;
			// const VType t31 = qy*t2*y;
			// const VType t32 = qx*t4*x;
			// const VType t33 = qx*t3*x;
			// const VType t34 = qy*t4*y;
			// const VType t35 = qy*t3*y;
			// const VType t36 = qw*qx*qz*y*2.0;
			// const VType t37 = 1.0/t6;
			// J( i*3+0, 0 ) = 1.0;
			// J( i*3+0, 1 ) = J( i*3+0, 2 ) = 0;
			// J( i*3+0, 3 ) = t7*(t11+t12+t13+t31+t34+t35+t36+qx*t2*x*2.0+qx*t3*x*2.0-qy*t5*y-qz*t5*z-qw*qx*qy*z*2.0)*2.0;
			// J( i*3+0, 4 ) = t7*(t8+t9+t10-qy*t4*x*2.0-qy*t5*x*2.0-qx*t2*y+qx*t3*y+qx*t4*y+qx*t5*y-qw*t2*z+qw*qy*qz*y*2.0-qx*qy*qz*z*2.0)*2.0;
			// J( i*3+0, 5 ) = t7*(t17+t18+t19+t21+t22+t24+qz*t4*x*2.0+qz*t5*x*2.0-qw*t3*y-qx*t2*z-qx*t4*z-qx*t5*z)*-2.0;
			// J( i*3+0, 6 ) = t7*(t14+t15+t16+t28+qw*t2*x*2.0+qw*t3*x*2.0-qz*t2*y-qz*t3*y-qz*t5*y-qy*t4*z-qw*qx*qy*y*2.0-qw*qx*qz*z*2.0)*2.0;
			// J( i*3+1, 1 ) = 1.0;
			// J( i*3+1, 0 ) = J( i*3+1, 2 ) = 0;
			// J( i*3+1, 3 ) = t7*(t8-t9+t10-qy*t2*x-qy*t3*x-qy*t4*x+qy*t5*x+qx*t2*y*2.0+qx*t4*y*2.0+qw*t2*z+qw*qx*qz*x*2.0+qx*qy*qz*z*2.0)*-2.0;
			// J( i*3+1, 4 ) = t7*(t11+t12-t13+t30+t32+t33-qx*t2*x+qy*t3*y*2.0+qy*t5*y*2.0+qz*t5*z-qw*qy*qz*x*2.0+qw*qx*qy*z*2.0)*2.0;
			// J( i*3+1, 5 ) = t7*(t14+t15-t16+t25+t26+t27+t29-qw*t3*x-qz*t2*y*2.0-qz*t4*y*2.0-qx*qy*qz*x*2.0+qw*qx*qz*z*2.0)*2.0;
			// J( i*3+1, 6 ) = t7*(t17+t18+t20+t23-qz*t2*x-qz*t3*x-qz*t5*x-qw*t3*y*2.0-qw*t5*y*2.0+qx*t2*z-qx*t4*z+qx*t5*z)*-2.0;
			// J( i*3+2, 2 ) = 1.0;
			// J( i*3+2, 0 ) = J( i*3+2, 1 ) = 0;
			// J( i*3+2, 3 ) = t7*(t19+t20-t21+t22+t23-t24+qz*t2*x+qz*t3*x-qz*t5*x+qw*t3*y-qx*t2*z*2.0-qx*t4*z)*2.0;
			// J( i*3+2, 4 ) = t7*(-t25-t26+t27+t28+t29-qw*t3*x-qz*t2*y+qz*t3*y+qz*t5*y+qy*t3*z*2.0+qy*t5*z*2.0-qx*qy*qz*x*2.0-qw*qx*qy*y*2.0)*2.0;
			// J( i*3+2, 5 ) = t7*(-t12+t30+t31+t32-t33+t34-t35-t36+qx*t2*x+qy*t5*y-qz*t2*z*2.0+qw*qy*qz*x*2.0)*2.0;
			// J( i*3+2, 6 ) = z*(qw*2.0+qw*t7*(-t2+t3+t5)*2.0)-qy*t37*x*2.0+qx*t37*y*2.0+qw*t7*x*(qw*qy*2.0-qx*qz*2.0)*2.0-qw*t7*y*(qw*qx*2.0+qy*qz*2.0)*2.0;
			
		}
		// std::cout << "input  :\n" << input << "\n";
		// std::cout << "J      :\n" << J.size1() << " x " << J.size2() << std::endl;
	}
};



template< typename InputIterator >
bool estimatePose6D_3D3D( const InputIterator iBeginA, const InputIterator iEndA
	, Math::Pose& pose
	, const InputIterator iBeginB, const InputIterator iEndB
	, const Math::Optimization::OptTerminate& criteria )
{
#ifndef HAVE_LAPACK
	return false;
#else
	typedef typename std::iterator_traits< InputIterator >::value_type vector_type;
	typedef typename vector_type::value_type T;
	
	
	const std::size_t n = std::distance( iBeginA, iEndA );
	assert( n > 2 );
	
	// commented 1 step: take input as initial pose, no estimation here
	// 1) estimate a first initial guess,
	// if( !estimatePose6D_3D3D( iBeginA, iEndA, pose, iBeginB, iEndB ) )
		// return false;
		
	// 2) prepare the expectation values of the minimization function
	Math::Vector< T > measurement = Math::Vector< T >::zeros( 3*n );
	std::size_t i = 0;
	for( InputIterator it( iBeginA ); it != iEndA; ++it, ++i )
	{
		measurement[ i*3+0 ] = (*it)[ 0 ];
		measurement[ i*3+1 ] = (*it)[ 1 ];
		measurement[ i*3+2 ] = (*it)[ 2 ];
	}
	
	// 3) set the evaluation function
	PoseEstimation3D3D::PointCorrespodencesSinglePose< InputIterator > func( iBeginB, iEndB );
	
	// 4) set the parameter vector to optimize
	// Math::Vector< T, 7 > paramVector;
	Math::Vector< T > paramVector( 7 );
	pose.toVector( paramVector );
	
	
	// 5) perform optimization
	T residual = Ubitrack::Math::Optimization::levenbergMarquardt( func, paramVector, measurement, criteria, Math::Optimization::OptNoNormalize() );	
	
	// pose = Math::Pose::fromVector( paramVector );
	// skip the upper version to normalize the pose directly
	
	pose = Math::Pose( 
		Math::Quaternion( paramVector[ 3 ], paramVector[ 4 ], paramVector[ 5 ], paramVector[ 6 ] ).normalize()
		, Math::Vector< T, 3 >( paramVector [ 0 ], paramVector [ 1 ], paramVector [ 2 ] )
		);
	
	return true;
#endif	// HAVE_LAPACK
}
	

UBITRACK_EXPORT bool estimatePose6D_3D3D( const std::vector< Math::Vector3f >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3f >& pointsB
	, const Math::Optimization::OptTerminate& criteria );
	

UBITRACK_EXPORT bool estimatePose6D_3D3D( const std::vector< Math::Vector3d >& pointsA
	, Math::Pose& pose
	, const std::vector< Math::Vector3d >& pointsB
	, const Math::Optimization::OptTerminate& criteria );


	
}}} // namespace Ubitrack::Algorithm::PoseEstimation3D3D

#endif //__UBITRACK_ALGROITHM_ABSOLUTE_ORIENTATION_OPTIMIZATION_H_INCLUDED__
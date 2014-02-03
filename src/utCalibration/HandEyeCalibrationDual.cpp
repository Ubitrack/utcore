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
 * Implementation of Dual Quaternion Hand-Eye Calibration
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */


#include "HandEyeCalibrationDual.h"

#include <utMath/Pose.h>
#include <utMath/Matrix.h>
#include <utMath/Blas1.h> // inner_product
#include <utMath/Stochastic/identity_iterator.h>

#include <algorithm> //std::transform

//shortcuts to namespaces
#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>
namespace lapack = boost::numeric::bindings::lapack;



namespace Ubitrack{ namespace Math {

template< typename T >
struct Pose2DualQuaternion
{
	typedef typename Math::Vector< T, 8 > return_value;
	
	Math::Vector< T, 8 > operator()( const Math::Pose& pose ) const
	{
		const T qw = pose.rotation().w();
		const T qx = pose.rotation().x();
		const T qy = pose.rotation().y();
		const T qz = pose.rotation().z();
		
		const T tx = pose.translation()[ 0 ];
		const T ty = pose.translation()[ 1 ];
		const T tz = pose.translation()[ 2 ];
		
		// dual-quaternion return value: (q | q')
		// order: q  :=  (qw  | qx  qy  qz )
		// order: q' :=  (q'w | q'x q'y q'z)
		Math::Vector< T, 8 > dualQuat;
		
		// first, the easy part: :)
		// quaternion goes into quaterion part (q) 
		dualQuat( 0 ) = qw;
		dualQuat( 1 ) = qx;
		dualQuat( 2 ) = qy;
		dualQuat( 3 ) = qz;
				
		// second, the tricky part:
		// translation t goes into quaterion dual part (q' == aka q prime )
		// how does it work?
		// you take translation t and make a quaternion q_t out of it,
		// assuming zero for the real part an (t/2) as the imaginary parts.
		// you apply quaternion multiplication to q_t with q from the right side
		// and done. 
		// summary: q' := q( t/2, 0) * q
		//
		// technically it is done in two parts
		// part 1: inner product of imaginary part from q and t
		// -(1/2) * (q_xyz' *t) 
		dualQuat( 4 ) = -0.5* ( tx*qx + ty*qy + tz*qz );
		
		// part 2:
		// (1/2) * (cross_product(t, q_xyz) + q_w*t)
		dualQuat( 5 ) = 0.5* ( ( ty*qz - tz*qy ) + qw*tx );
		dualQuat( 6 ) = 0.5* ( ( tz*qx - tx*qz ) + qw*ty );
		dualQuat( 7 ) = 0.5* ( ( tx*qy - ty*qx ) + qw*tz );
		
		return dualQuat;
	}
};

template< typename T >
struct difference_dual_a
{
	Math::Vector< T, 8 > operator()( const Math::Pose& pose1, const Math::Pose& pose2 ) const
	{
		const Math::Pose pose = (~pose2) * pose1;
		return Pose2DualQuaternion< T >()( pose );
	}
};

template< typename T >
struct difference_dual_b
{		
	Math::Vector< T, 8 > operator()( const Math::Pose& pose1, const Math::Pose& pose2 ) const
	{
		const Math::Pose pose = (pose2) * (~pose1);
		return Pose2DualQuaternion< T >()( pose );
	}
};

struct SolveQuadratic
{
	template< typename T >
	Math::Vector< T, 2 > operator()( const Math::Vector< T, 3 >& quadratic ) const
	{
		// solve ax^2 + bx + c = 0 for x_{1,2}
		// assume input vector with following order[ c b a ]
		const T c = quadratic[ 0 ];
		const T b = quadratic[ 1 ];
		const T m2a = 2.*quadratic[ 2 ]; // <-- 2a is needed anyway
		
		const T root = std::sqrt( (b*b)-(2.*m2a*c) );// <- 2a instead of 4a, see above
		const T x1 = (-b + root) / m2a;// <- here as well 2a instead of 4a
		const T x2 = -(b + root) / m2a;
		return Math::Vector< T, 2 >( x1, x2 );
	}
};

}} //namesapace Ubitrack::Math;

namespace Ubitrack { namespace Calibration {

namespace {


/** 
 * @internal
 * @brief 
 * This function computes a result among the input values from first to last using the given binary operator
 *
 * This template function is inspired from http://www.cplusplus.com/reference/numeric/adjacent_difference/
 */
template <class InputIterator, class OutputIterator, class BinaryOperation >
OutputIterator own_adjacent_difference ( InputIterator first, InputIterator last, OutputIterator result, BinaryOperation binary_op )
{
	if (first!=last)
	{
		typename std::iterator_traits< InputIterator >::value_type val,prev;
		prev = *first;
		while(++first!=last)
		{
			val = *first;
			*++result++ = binary_op( val,prev );
			prev = val;
		}
		++result;
	}
	return result;
}

template< typename ForwardIterator >
void generate_diff_poses( const ForwardIterator begin, const ForwardIterator end )
{

}

template< bool use_all_pairs, typename T, typename ForwardIterator >
bool estimatePose6D_6D6D_impl( const ForwardIterator itBeginEye, const ForwardIterator itEndEye,
	Math::Pose& pose,
	const ForwardIterator itBeginHand, const ForwardIterator itEndHand )
{
	typedef typename Ubitrack::Math::Pose2DualQuaternion< T >::return_value dual_type;
	
	const std::size_t n_in = std::distance( itBeginEye, itEndEye );
	assert( n_in > 2 );
	assert( n_in == std::distance( itBeginHand, itEndHand ) );
	const std::size_t n = n_in-1;
	const std::size_t m = use_all_pairs ? (n*n_in)/2 : n;
	
	// std::cout << "Received " << n_in << " poses, using " << m << " difference poses.\n";
	
	std::vector< dual_type > dualA; // <- paper from Daniilidis uses a and b to specify dual quaternions
	{ // set the first dual quaternion from input differences
		dualA.reserve( m );
		if( !use_all_pairs )
			own_adjacent_difference ( itBeginEye, itEndEye, std::back_inserter( dualA ), Ubitrack::Math::difference_dual_a< T >() );
		else
		{
			ForwardIterator itPose = itBeginEye;
			for( ;itPose != itEndEye; )
			{
				Ubitrack::Util::identity< Math::Pose > id( *itPose );
				std::advance( itPose, 1 );
				std::transform( itPose , itEndEye, id.begin(), std::back_inserter( dualA ), Ubitrack::Math::difference_dual_a< T >() );
			}
		}
	}
	
	std::vector< dual_type > dualB;  // <- paper from Daniilidis uses a and b to specify dual quaternions
	{
		dualB.reserve( m );
		if( !use_all_pairs )
			own_adjacent_difference ( itBeginHand, itEndHand, std::back_inserter( dualB ), Ubitrack::Math::difference_dual_b< T >() );
		else
		{
			ForwardIterator itPose = itBeginHand;
			for( ;itPose != itEndHand; )
			{
				Ubitrack::Util::identity< Math::Pose > id( *itPose );
				std::advance( itPose, 1 );
				std::transform( itPose , itEndHand, id.begin(), std::back_inserter( dualB ), Ubitrack::Math::difference_dual_b< T >() );
			}
		}
	}
	
	assert( dualA.size() == dualB.size() );
	assert( dualA.size() == m );
	
	// create a 6*n-by-8 matrix with the following scheme:
	// [a  - b ] [a  + b ]_x [0 0 0]^T [0_{3x3}]
	// [a' - b'] [a' + b']_x [a - b]   [a + b ]_x
	Math::Matrix< T > matrixT( 6*m, 8 ); // <- all values get in there
	for( std::size_t index = 0; index<m; ++index )
	{
		const dual_type a = dualA[ index ];
		const dual_type b = dualB[ index ];

		const T ax = a[ 1 ];
		const T ay = a[ 2 ];
		const T az = a[ 3 ];
		const T apx = a[ 5 ];
		const T apy = a[ 6 ];
		const T apz = a[ 7 ];
		
		const T bx = b[ 1 ];
		const T by = b[ 2 ];
		const T bz = b[ 3 ];
		const T bpx = b[ 5 ];
		const T bpy = b[ 6 ];
		const T bpz = b[ 7 ];

		const std::size_t i = (index*6);
		matrixT( i+0, 0 ) = ax-bx;
		matrixT( i+1, 0 ) = ay-by;
		matrixT( i+2, 0 ) = az-bz;
		matrixT( i+3, 0 ) = apx-bpx;
		matrixT( i+4, 0 ) = apy-bpy;
		matrixT( i+5, 0 ) = apz-bpz;
		
		// second column (first column cross product matrix)
		matrixT( i+0, 1 ) = 0;
		matrixT( i+1, 1 ) = az+bz;
		matrixT( i+2, 1 ) = -(ay+by);
		matrixT( i+3, 1 ) = 0;
		matrixT( i+4, 1 ) = apz+bpz;
		matrixT( i+5, 1 ) = -(apy+bpy);
		
		matrixT( i+0, 2 ) = -(az+bz);
		matrixT( i+1, 2 ) = 0;
		matrixT( i+2, 2 ) = ax+bx;
		matrixT( i+3, 2 ) = -(apz+bpz);
		matrixT( i+4, 2 ) = 0 ;
		matrixT( i+5, 2 ) = apx+bpx ;
		
		matrixT( i+0, 3 ) = (ay+by);
		matrixT( i+1, 3 ) = -(ax+bx);
		matrixT( i+2, 3 ) = 0;
		matrixT( i+3, 3 ) = (apy+bpy);
		matrixT( i+4, 3 ) = -(apx+bpx);
		matrixT( i+5, 3 ) = 0;
		
		matrixT( i+0, 4 ) = matrixT( i+1, 4 ) = matrixT( i+2, 4 ) = 0;
		matrixT( i+3, 4 ) = ax-bx;
		matrixT( i+4, 4 ) = ay-by;
		matrixT( i+5, 4 ) = az-bz;
		
		matrixT( i+0, 5 ) = matrixT( i+1, 5 ) = matrixT( i+2, 5 ) = 0;
		matrixT( i+3, 5 ) = 0;
		matrixT( i+4, 5 ) = az+bz;
		matrixT( i+5, 5 ) = -(ay+by);
		
		matrixT( i+0, 6 ) = matrixT( i+1, 6 ) = matrixT( i+2, 6 ) = 0;
		matrixT( i+3, 6 ) = -(az+bz);
		matrixT( i+4, 6 ) = 0;
		matrixT( i+5, 6 ) = ax+bx;
		 
		matrixT( i+0, 7 ) = matrixT( i+1, 7 ) = matrixT( i+2, 7 ) = 0;
		matrixT( i+3, 7 ) = (ay+by);
		matrixT( i+4, 7 ) = -(ax+bx);
		matrixT( i+5, 7 ) = 0;
	}
		
	Math::Vector< T, 8 > s;
	Math::Matrix< T > u( 6*m, 6*m );
	Math::Matrix< T, 8, 8 > vt;
	
	// see http://dlib.net/dlib/matrix/lapack/gesvd.h.html for parameters of 'gesvd'
	int i = boost::numeric::bindings::lapack::gesvd( 'N', 'S', matrixT, s, u, vt );
	if( i != 0 )
	{
		std::cout << "SVD returned with error " << i << "\n";
		return false;
	}
	
	
	// check if last two singular values are smallest (near zero)
	// and other ones are bigger. 
	// luckily lapack returns smallest values always at the end
	const T epsilon( 1e-02 );
	if( s( 7 ) > epsilon || s( 6 ) > epsilon || s( 5 ) < epsilon )
	{
		std::cout << "Matrix Vt\n" << vt << "\n and corresponding singular values.\n";
		std::cout << "Check the singular values for debugging:\n" << s << "\n";	
		return false;
	}

	Math::Vector< T, 4 > u1( vt( 6, 0 ), vt( 6, 1 ), vt( 6, 2 ), vt( 6, 3 ) );
	Math::Vector< T, 4 > v1( vt( 6, 4 ), vt( 6, 5 ), vt( 6, 6 ), vt( 6, 7 ) );
	Math::Vector< T, 4 > u2( vt( 7, 0 ), vt( 7, 1 ), vt( 7, 2 ), vt( 7, 3 ) );
	Math::Vector< T, 4 > v2( vt( 7, 4 ), vt( 7, 5 ), vt( 7, 6 ), vt( 7, 7 ) );
	
	const T a = inner_product( u1, v1 );
	const T b = inner_product( u1, v2 ) + inner_product( u2, v1 ); 
	const T c = inner_product( u2, v2 );
	const Math::Vector< T, 3 > quEq( c, b, a );
	const Math::Vector< T, 2 > s12 = Ubitrack::Math::SolveQuadratic()( quEq );
	
	const T dotU1 = inner_product( u1, u1 );
	const T dotU1U2_2 = inner_product( u1, u2 ) * 2; // <- 2 times inner product
	const T dotU2 = inner_product( u2, u2 );
	const T s1 = s12[ 0 ] * s12[ 0 ] * dotU1 + s12[ 0 ] * dotU1U2_2 + dotU2;
	const T s2 = s12[ 1 ] * s12[ 1 ] * dotU1 + s12[ 1 ] * dotU1U2_2 + dotU2;
	const T lambda2 = (s1 > s2) ? std::sqrt( 1/s1 ) : std::sqrt( 1/s2 );
	const T lambda1 = (s1 > s2) ? lambda2 * s12[ 0 ] : lambda2 * s12[ 1 ];
	
	if( (lambda1 != lambda1) || (lambda2 != lambda2) )
	{
		std::cout << "No valid parameters for lambda 1/2 = " << lambda1 << " " << lambda2 << "\n";
		return false;
	}

	// prepare the result
	Math::Vector< T, 4 > q = (lambda1 * u1) + (lambda2 * u2); 
	Math::Vector< T, 4 > qp = (lambda1 * v1) + (lambda2 * v2);
	Math::Quaternion qprime( qp( 1 ), qp( 2 ), qp( 3 ), qp( 0 ) );
	Math::Quaternion q_conj( -q( 1 ), -q( 2 ), -q( 3 ), q( 0 ) );
	Math::Quaternion t_final = qprime*q_conj;
	
	pose = Ubitrack::Math::Pose( Ubitrack::Math::Quaternion( q( 1 ), q( 2 ), q( 3 ), q( 0 )), Ubitrack::Math::Vector< T, 3 >( 2*t_final.x(), 2*t_final.y(), 2*t_final.z() ) );
	
	return true;
};
	
} // anonymous-namespace

UBITRACK_EXPORT bool estimatePose6D_6D6D( const std::vector< Math::Pose >& eyes, Math::Pose& pose,
	const std::vector< Math::Pose >& hands )
{
	// // uses all distinct pairs of input poses:
	return estimatePose6D_6D6D_impl< true, double >( eyes.begin(), eyes.end(), pose, hands.begin(), hands.end() );
	
	// uses only the direct paris from input poses:
	// return estimatePose6D_6D6D_impl< false, double >( eyes.begin(), eyes.end(), pose, hands.begin(), hands.end() );
}



}} // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

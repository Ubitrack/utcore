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


#include "DualQuaternion.h"
#include "DataSelection.h"

#include <utMath/Matrix.h>
#include <utMath/Blas1.h> // inner_product

#include <algorithm> //std::transform

//shortcuts to namespaces
#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>

namespace Ubitrack{ namespace Math {

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

namespace Ubitrack { namespace Algorithm { namespace PoseEstimation6D6D {

namespace {

template< typename InputIterator >
bool estimatePose6D_6D6D_impl( InputIterator itBeginEye, const InputIterator itEndEye,
	Math::Pose& pose,
	InputIterator itBeginHand, const InputIterator itEndHand )
{
	typedef typename std::iterator_traits< InputIterator >::value_type dual_type;
	typedef typename dual_type::value_type value_type;
	
	const std::size_t n = std::distance( itBeginEye, itEndEye );
	const std::size_t n2 = std::distance( itBeginHand, itEndHand );
	assert( n == n2 );
	assert( n > 2 );
	
	// create a 6*n-by-8 matrix with the following scheme:
	// [a  - b ] [a  + b ]_x [0 0 0]^T [0_{3x3}]
	// [a' - b'] [a' + b']_x [a - b]   [a + b ]_x
	Math::Matrix< value_type > matrixT( 6*n, 8 ); // <- all values get in there
	for( std::size_t index = 0; index<n; ++index, ++itBeginEye, ++itBeginHand )
	{
		const dual_type a = *itBeginEye;
		const dual_type b = *itBeginHand;

		const value_type ax = a[ 1 ];
		const value_type ay = a[ 2 ];
		const value_type az = a[ 3 ];
		const value_type apx = a[ 5 ];
		const value_type apy = a[ 6 ];
		const value_type apz = a[ 7 ];
		
		const value_type bx = b[ 1 ];
		const value_type by = b[ 2 ];
		const value_type bz = b[ 3 ];
		const value_type bpx = b[ 5 ];
		const value_type bpy = b[ 6 ];
		const value_type bpz = b[ 7 ];

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
		
	Math::Vector< value_type, 8 > s;
	Math::Matrix< value_type > u( 6*n, 6*n );
	Math::Matrix< value_type, 8, 8 > vt;
	
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
	const value_type epsilon( 1e-02 );
	if( s( 7 ) > epsilon || s( 6 ) > epsilon || s( 5 ) < epsilon )
	{
		std::cout << "Matrix Vt\n" << vt << "\n and corresponding singular values.\n";
		std::cout << "Check the singular values for debugging:\n" << s << "\n";	
		return false;
	}

	Math::Vector< value_type, 4 > u1( vt( 6, 0 ), vt( 6, 1 ), vt( 6, 2 ), vt( 6, 3 ) );
	Math::Vector< value_type, 4 > v1( vt( 6, 4 ), vt( 6, 5 ), vt( 6, 6 ), vt( 6, 7 ) );
	Math::Vector< value_type, 4 > u2( vt( 7, 0 ), vt( 7, 1 ), vt( 7, 2 ), vt( 7, 3 ) );
	Math::Vector< value_type, 4 > v2( vt( 7, 4 ), vt( 7, 5 ), vt( 7, 6 ), vt( 7, 7 ) );
	
	const value_type a = inner_product( u1, v1 );
	const value_type b = inner_product( u1, v2 ) + inner_product( u2, v1 ); 
	const value_type c = inner_product( u2, v2 );
	const Math::Vector< value_type, 3 > quEq( c, b, a );
	const Math::Vector< value_type, 2 > s12 = Ubitrack::Math::SolveQuadratic()( quEq );
	
	const value_type dotU1 = inner_product( u1, u1 );
	const value_type dotU1U2_2 = inner_product( u1, u2 ) * 2; // <- 2 times inner product
	const value_type dotU2 = inner_product( u2, u2 );
	const value_type s1 = s12[ 0 ] * s12[ 0 ] * dotU1 + s12[ 0 ] * dotU1U2_2 + dotU2;
	const value_type s2 = s12[ 1 ] * s12[ 1 ] * dotU1 + s12[ 1 ] * dotU1U2_2 + dotU2;
	const value_type lambda2 = (s1 > s2) ? std::sqrt( 1/s1 ) : std::sqrt( 1/s2 );
	const value_type lambda1 = (s1 > s2) ? lambda2 * s12[ 0 ] : lambda2 * s12[ 1 ];
	
	if( (lambda1 != lambda1) || (lambda2 != lambda2) )
	{
		std::cout << "No valid parameters for lambda 1/2 = " << lambda1 << " " << lambda2 << "\n";
		return false;
	}

	// prepare the result
	Math::Vector< value_type, 4 > q = (lambda1 * u1) + (lambda2 * u2); 
	Math::Vector< value_type, 4 > qp = (lambda1 * v1) + (lambda2 * v2);
	Math::Quaternion qprime( qp( 1 ), qp( 2 ), qp( 3 ), qp( 0 ) );
	Math::Quaternion q_conj( -q( 1 ), -q( 2 ), -q( 3 ), q( 0 ) );
	Math::Quaternion t_final = qprime*q_conj;
	
	pose = Ubitrack::Math::Pose( Ubitrack::Math::Quaternion( q( 1 ), q( 2 ), q( 3 ), q( 0 )), Ubitrack::Math::Vector< value_type, 3 >( 2*t_final.x(), 2*t_final.y(), 2*t_final.z() ) );
	
	return true;
};
	
} // anonymous-namespace

bool estimatePose6D_6D6D( const std::vector< Math::Pose >& eyes, Math::Pose& pose,
	const std::vector< Math::Pose >& hands )
{

	typedef Ubitrack::Math::Vector< double, 8 > dual_type; // defining the input values of the algorithm.
	const bool use_all_pairs = true; // <-- sign if all distinct pairs of poses should be used or only the neighbouring

	
	// checking the validity of inputs
	const std::size_t n_in = std::distance( eyes.begin(), eyes.end() );
	const std::size_t n_in2 = std::distance( hands.begin(), hands.end() );
	assert( n_in == n_in2 );
	assert( n_in > 2 ); // <- algorithm needs at least 3 relative movements
	
	// determining the amount of relative pose movements
	const std::size_t n = n_in-1;
	const std::size_t m = use_all_pairs ? (n*n_in)/2 : n;
	
	
	// generate the relative pose movements for the eye in forward direction
	std::vector< dual_type > dualA; // <- paper from Daniilidis uses a and b to specify dual quaternions
	dualA.reserve( m );
	generate_relative_pose6D_impl< use_all_pairs, true >( eyes.begin(), eyes.end(), std::back_inserter( dualA ) );
	
	// generate the relative pose movements for the hand in backward direction
	std::vector< dual_type > dualB;  // <- paper from Daniilidis uses a and b to specify dual quaternions
	dualB.reserve( m );
	generate_relative_pose6D_impl< use_all_pairs, false >( hands.begin(), hands.end(), std::back_inserter( dualB ) );

	return estimatePose6D_6D6D_impl( dualA.begin(), dualA.end(), pose, dualB.begin(), dualB.end() );
	
}

/// @internal overloaded function that takes dual quaternions (as 8-vector) assuming them to be relative poses.
bool estimatePose6D_6D6D( const std::vector< Math::Vector< double, 8 > >& eyes, Math::Pose& pose,
	const std::vector< Math::Vector< double, 8 > >& hands )
{
	return estimatePose6D_6D6D_impl( eyes.begin(), eyes.end(), pose, hands.begin(), hands.end() );
}

}}} // namespace Ubitrack::Algorithm::PoseEstimation6D6D

#endif // HAVE_LAPACK


#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointTransformation.h>
#include <utAlgorithm/PoseEstimation3D3D/AbsoluteOrientation.h>


#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;


void fillDemoVectorsDeterministic( Vector< double, 3 >* left, Vector< double, 3 >* right, Quaternion q, Vector< double, 3 > t )
{
	left[0][0] = 1.0;
	left[0][1] = 0.0;
	left[0][2] = 0.0;

	left[1][0] = 0.0;
	left[1][1] = 1.0;
	left[1][2] = 0.0;

	left[2][0] = 0.0;
	left[2][1] = 0.0;
	left[2][2] = 1.0;

	left[3][0] = 0.5;
	left[3][1] = 0.5;
	left[3][2] = 0.7;

	right[0] = q*left[0]+t;
	right[1] = q*left[1]+t;
	right[2] = q*left[2]+t;
	right[3] = q*left[3]+t;
}

void testAbsoluteOrientationDeterministic()
{
	// Vector< double, 3 > axis ( 1.0, 1.0, 1.5 );
	// Quaternion q ( axis, M_PI/6.0 );
	// Vector< double, 3 > t ( -1.0, 3.0, 2.5 );

	// Vector< double, 3 > left[4];
	// Vector< double, 3 > right[4];
	// fillDemoVectorsDeterministic ( left, right, q, t );

	// Pose p = Ubitrack::Algorithm::calculateAbsoluteOrientation ( &left[0], &left[4], &right[0], &right[4] );

	// BOOST_CHECK_SMALL ( vectorDiff ( p.translation(), t ), 10e-6 );
	// BOOST_CHECK_SMALL ( quaternionDiff ( p.rotation(), q ), 10e-6 );
}

template< typename T >
void testAbsoluteOrientationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -10, 10 );
	// typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.01, 0.01 ); // uniform noise for translation
	typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.005 ); // gaussian noise for translation
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n_p3d = 3 + ( iRun % 28);//( Random::distribute_uniform< std::size_t >( 3, 30 ) );

		std::vector< Vector< T, 3 > > rightFrame;
		rightFrame.reserve( n_p3d );
		std::generate_n ( std::back_inserter( rightFrame ), n_p3d,  randVector );
		
		
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		Matrix< T, 3, 4 > trafo( q, t );
		
		std::vector< Vector< T, 3 > > leftFrame;
		leftFrame.reserve( n_p3d );
		Geometry::transform_points( trafo, rightFrame.begin(), rightFrame.end(), std::back_inserter( leftFrame ) );

		// do some estimation now
		Pose estimatedPose;
		const bool b_done = Ubitrack::Algorithm::PoseEstimation3D3D::estimatePose6D_3D3D( leftFrame, estimatedPose, rightFrame );

		BOOST_WARN_MESSAGE( b_done, "Algorithm did not successfully estimate a result with " << n_p3d << " points." );
		
		
		if( !b_done )
			continue;
	
		{	// calculate some errors on resulting transformation
			const T rotDiff = quaternionDiff( estimatedPose.rotation(), q );
			const T posDiff = vectorDiff( estimatedPose.translation(), t );
			
			// check if pose is better than before (only for valid results)
			BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nCompare translation estimation using " << n_p3d << " points, error=" << posDiff  << ":\n" << t << " (expected)\n" << estimatedPose.translation() << " (estimated)\n");
			BOOST_CHECK_MESSAGE( rotDiff < epsilon, "\nCompare rotation estimation using " << n_p3d << " points, error=" << rotDiff  << ":\n" << q << " (expected)\n" << estimatedPose.rotation() << " (estimated)\n" );
		}
		
		{	// calculate an residual error from input data
			// const T err = Ubitrack::Algorithm::PoseEstimation3D3D::estimatePose6DResidual< T >( leftFrame.begin(), leftFrame.end(), estimatedPose, rightFrame.begin(), rightFrame.end() );
			// std::cout << "Residual Error is " << err.first << " stdDev: " << err.second << std::endl;
			// BOOST_CHECK_MESSAGE( err.first < epsilon, "\nResidual error using " << n_p3d << " points is too high :\n" << epsilon << " (expected)\n" << err.first << " (estimated)\n");
		}
	}
	
}
#ifndef HAVE_LAPACK
void TestAbsoluteOrientation()
{
	// Absolute Orientation does not work without lapack
}
#else // HAVE_LAPACK

void TestAbsoluteOrientation()
{
	// first do a deterministic test
	// testAbsoluteOrientationDeterministic();
	
	// do some iterations of random tests
	testAbsoluteOrientationRandom< float >( 10000, 1e-2f );
	testAbsoluteOrientationRandom< double >( 10000, 1e-6 );
}

#endif // HAVE_LAPACK

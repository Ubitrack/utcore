
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointTransformation.h>
#include <utAlgorithm/PoseEstimation3D3D/Ransac.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

template< typename T >
void testRansacAbsoluteOrientationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -100, 100 );
	// typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.01, 0.01 ); // uniform noise for translation
	// typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.005 ); // gaussian noise for translation
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n_p3d = 10+(iRun%91);

		std::vector< Vector< T, 3 > > rightFrame;
		rightFrame.reserve( n_p3d );
		std::generate_n ( std::back_inserter( rightFrame ), n_p3d,  randVector );
		
		
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		Matrix< T, 3, 4 > trafo( q, t );
		
		std::vector< Vector< T, 3 > > leftFrame;
		leftFrame.reserve( n_p3d );
		Geometry::transform_points( trafo, rightFrame.begin(), rightFrame.end(), std::back_inserter( leftFrame ) );
		
		// {	// add noise to the points (maybe only in optimized version
			// // add some noise to left side points
			// std::vector< Vector< T, 3 > > points_3D_noisy_left;
			// points_3D_noisy_left.reserve( n_p3d );
			// std::generate_n ( std::back_inserter( points_3D_noisy_left ), n_p3d,  randPositionNoise );
			// std::transform( leftFrame.begin(), leftFrame.end(), points_3D_noisy_left.begin(), leftFrame.begin(), std::plus< Vector< T, 3 > >() );
			
			// // add some noise to right side points
			// std::vector< Vector< T, 3 > > points_3D_noisy_right;
			// points_3D_noisy_right.reserve( n_p3d );
			// std::generate_n ( std::back_inserter( points_3D_noisy_right ), n_p3d,  randPositionNoise );
			// std::transform( rightFrame.begin(), rightFrame.end(), points_3D_noisy_right.begin(), rightFrame.begin(), std::plus< Vector< T, 3 > >() );
		// }
		
		
		// now produce some (10%) outlier on both sides
		const std::size_t outlier( n_p3d/10 );
		for( std::size_t i = 0; i<outlier; ++i )
		{
			const std::size_t index1 = Random::distribute_uniform< std::size_t >( 0, n_p3d-1 ) ;
			leftFrame[ index1 ] = randVector();
			
			const std::size_t index2 = Random::distribute_uniform< std::size_t >( 0, n_p3d-1 ) ;
			rightFrame[ index2 ] = randVector();
		}
		
		// set reasonable? ransac parameters
		const T threshold( 0.05 ); // <- this works quite well here, since there are some big outlier, might be different in other scenarios
		const std::size_t minSetSize( 3 ); // <- the algorithm needs at least 3 poses to estimate a result (minimum consensus set)
		const T percentOutlier( 0.4 ); // <- should be a reasonable number, maybe derived empirically
		const T sucessProbability( 0.99 ); // <- usually you want to have this parameter quite high (~95 or higher)
		const Ubitrack::Math::Optimization::RansacParameter< T > params( threshold, minSetSize, n_p3d, percentOutlier, sucessProbability );
		
		//BOOST_MESSAGE( "Starting with at least " << params.nMaxIterations << " iterations for " << params.nMinInlier << " expected inlier of " << n << " poses."  );
		
		// do some robust estimation now
		Pose estimatedPose;
		const bool b_done = Ubitrack::Algorithm::PoseEstimation3D3D::estimatePose6D_3D3D( leftFrame, estimatedPose, rightFrame, params );
		
		// calculate some errors
		const T rotDiff = quaternionDiff( estimatedPose.rotation(), q );
		const T posDiff = vectorDiff( estimatedPose.translation(), t );
		
		if( !b_done )
		{
			BOOST_WARN_MESSAGE( b_done, "Algorithm did not successfully estimate a result with " << n_p3d 
				<< " points.\nRemaining difference in rotation " << rotDiff << ", difference in translation " << posDiff << "." );
			continue;
		}
		
		// check if pose is better than before (only for valid results)
		BOOST_CHECK_MESSAGE( rotDiff < epsilon, "\nCompare rotation    result (expected vs. estimated) using " << n_p3d << " points:\n" << q << " " << estimatedPose.rotation() );
		BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nCompare translation result (expected vs. estimated) using " << n_p3d << " points:\n" << t << " " << estimatedPose.translation() );
	}	
}

#ifndef HAVE_LAPACK
void TestRobustAbsoluteOrientation()
{
	// Absolute Orientation does not work without lapack
}

#else // HAVE_LAPACK

void TestRobustAbsoluteOrientation()
{
	// do some iterations of random tests
	testRansacAbsoluteOrientationRandom< float >( 10000, 1e-2f );
	testRansacAbsoluteOrientationRandom< double >( 10000, 1e-6 );
}

#endif // HAVE_LAPACK

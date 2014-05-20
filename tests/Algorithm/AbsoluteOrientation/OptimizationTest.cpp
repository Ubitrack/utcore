
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointTransformation.h>
#include <utAlgorithm/PoseEstimation3D3D/AbsoluteOrientation.h>
#include <utAlgorithm/PoseEstimation3D3D/Optimization.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

template< typename T >
void testOptimizedAbsoluteOrientationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -1, 1 );
	typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.01, 0.01 ); // uniform noise for translation
	// typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.01 ); // gaussian noise for translation
	
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
		
		// perturb the pose little bit
		const T rotEps( 0.015 );
		const Pose perurbedPose( Quaternion( q.x() + Random::distribute_uniform< T >( -rotEps, rotEps )
				, q.y() + Random::distribute_uniform< T >( -rotEps, rotEps )
				, q.z() + Random::distribute_uniform< T >( -rotEps, rotEps )
				, q.w() + Random::distribute_uniform< T >( -rotEps, rotEps ) ).normalize()
				, t + randPositionNoise() );
		
		const T rotDiff1 = quaternionDiff( perurbedPose.rotation(), q );
		const T posDiff1 = vectorDiff( perurbedPose.translation(), t );
		
		// set reasonable? optimization parameters
		const Ubitrack::Math::Optimization::OptTerminate termCrit( 50, 1e-8 );	
		
		Pose optimizedPose = perurbedPose;
		const bool b_done2 = Ubitrack::Algorithm::PoseEstimation3D3D::estimatePose6D_3D3D( leftFrame, optimizedPose, rightFrame, termCrit );
		
		// calculate some errors
		const T rotDiff2 = quaternionDiff( optimizedPose.rotation(), q );
		const T posDiff2 = vectorDiff( optimizedPose.translation(), t );
		
		if( !b_done2 )
		{
			BOOST_WARN_MESSAGE( b_done2, "Algorithm did not successfully estimate a result with " << n_p3d 
				<< " points.\nRemaining difference in rotation " << rotDiff2 << ", difference in translation " << posDiff2 << "." );
			continue;
		}
		
		// check if optimized pose is better than other one
		BOOST_CHECK_MESSAGE( rotDiff2 < rotDiff1, "\nCompare rotation    result using " << n_p3d << " points:\n" << q << " (expected)\n" << perurbedPose.rotation() << " (perturbed)\n" << optimizedPose.rotation() << " (optimized)" );
		BOOST_CHECK_MESSAGE( posDiff2 < posDiff1, "\nCompare translation result using " << n_p3d << " points:\n" << t << " (expected)\n" << perurbedPose.translation() << " (perturbed)\n" << optimizedPose.translation() << " (optimized)" );
	}	
}

#ifdef HAVE_LAPACK

void TestOptimizedAbsoluteOrientation()
{
	// do some iterations of random tests
	testOptimizedAbsoluteOrientationRandom< float >( 1000, 1e-2f );
	testOptimizedAbsoluteOrientationRandom< double >( 1000, 1e-6 );
}

#else // HAVE_LAPACK

void TestRobustAbsoluteOrientation()
{
	// Absolute Orientation does not work without lapack
}

#endif // HAVE_LAPACK

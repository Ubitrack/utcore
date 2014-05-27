
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointTransformation.h>
#include <utMath/Stochastic/Average.h>
#include <utMath/Stochastic/Gaussian.h>
#include <utAlgorithm/PoseEstimation3D3D/AbsoluteOrientation.h>
#include <utAlgorithm/PoseEstimation3D3D/ErrorEstimation.h>
#include <utAlgorithm/PoseEstimation3D3D/CovarianceEstimation.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

template< typename T >
void testCovarianceAbsoluteOrientationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -10, 10 );
	// typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.01, 0.01 ); // uniform noise for translation
	typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.05 ); // gaussian noise for translation
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n_p3d = 3 + ( iRun % 28);//( Random::distribute_uniform< std::size_t >( 3, 30 ) );

		std::vector< Vector< T, 3 > > rightFrame;
		rightFrame.reserve( n_p3d );
		std::generate_n ( std::back_inserter( rightFrame ), n_p3d,  randVector );
		
		
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		Pose pose ( q, t );
		Matrix< T, 3, 4 > trafo( q, t );
		
		// stores left hand results
		std::vector< Vector< T, 3 > > leftFrame;
		leftFrame.reserve( n_p3d );
		
		// stores noisy pose data for perturbation
		std::vector< Pose > noisyPoses;
		noisyPoses.reserve( n_p3d );
		
		{	// generate some random poses via random rotations
			
			for( std::size_t i = 0; i<n_p3d; ++i )
			{
				const T rotEps( 0.15 );
				const Pose noisyPose( Quaternion( q.x() + Random::distribute_uniform< T >( -rotEps, rotEps )
					, q.y() + Random::distribute_uniform< T >( -rotEps, rotEps )
					, q.z() + Random::distribute_uniform< T >( -rotEps, rotEps )
					, q.w() + Random::distribute_uniform< T >( -rotEps, rotEps ) ).normalize()
					, t + randPositionNoise() );
				
				noisyPoses.push_back( noisyPose );
				leftFrame.push_back( noisyPose * rightFrame[ i ] );
			}
		}
		// Stochastic::Gaussian< T, 7 > gauss;
		Ubitrack::Math::Stochastic::Average< Ubitrack::Math::Pose, Ubitrack::Math::ErrorPose > avg;
		Ubitrack::Math::ErrorPose errPose = avg.mean( noisyPoses );
		std::cout << "Error Pose " << errPose << "\n";
		// bool estimate_gaussian( const InputIterator itBegin, const InputIterator itEnd, Stochastic::Gaussian< T, N > &gaussian )
		
		// {	// add noise to the points (only if you like to)
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
		
		
		// do some estimation now
		Pose estimatedPose;
		const bool b_done = Ubitrack::Algorithm::PoseEstimation3D3D::estimatePose6D_3D3D( leftFrame, estimatedPose, rightFrame );

		BOOST_WARN_MESSAGE( b_done, "Algorithm did not successfully estimate a result with " << n_p3d << " points." );
		
		
		if( !b_done )
			continue;
	
		
		{	// calculate an residual error from input data
			// const T err = Ubitrack::Algorithm::PoseEstimation3D3D::estimatePose6DResidual< T >( leftFrame.begin(), leftFrame.end(), estimatedPose, rightFrame.begin(), rightFrame.end() );
			// std::cout << "Residual Error is " << err.first << " stdDev: " << err.second << std::endl;
			// BOOST_CHECK_MESSAGE( err.first < epsilon, "\nResidual error using " << n_p3d << " points is too high :\n" << epsilon << " (expected)\n" << err.first << " (estimated)\n");
		}
		
		{	// estimate a covariance to the corresponding estimated pose
			Ubitrack::Math::Matrix< T, 6, 6 > cov;
			const bool res = Ubitrack::Algorithm::PoseEstimation3D3D::estimatePose6DCovariance( leftFrame.begin(), leftFrame.end(), estimatedPose, rightFrame.begin(), rightFrame.end(), cov );
			
			const T posErr = std::sqrt( cov( 0, 0 ) + cov( 1, 1 ) +  cov( 2, 2 ) );
			const T rotErr = std::sqrt( cov( 3, 3 ) + cov( 4, 4 ) +  cov( 5, 5 ) );
			std::cout << "Pose Covariance " << cov << "\n";
			// std::cout << "Translation Error using " << n_p3d << " points:\n" << std::sqrt( cov( 0, 0 ) + cov( 1, 1 ) +  cov( 2, 2 ) )<< "\n";
			// std::cout << "Quaternion  Error using " << n_p3d << " points:\n" << std::sqrt( cov( 3, 3 ) + cov( 4, 4 ) +  cov( 5, 5 ) )<< "\n";
			// BOOST_CHECK_MESSAGE( posErr < epsilon, "\nCompare translation estimation using " << n_p3d << " points, stdDev=" << posErr  << ":\n" << t << " (expected)\n" << estimatedPose.translation() << " (estimated)\n");
			// BOOST_CHECK_MESSAGE( rotErr < epsilon, "\nCompare rotation estimation using " << n_p3d << " points, stdDev=" << rotErr  << ":\n" << q << " (expected)\n" << estimatedPose.rotation() << " (estimated)\n" );
		}
	}
	
}
#ifndef HAVE_LAPACK
void TestCovarianceAbsoluteOrientation()
{
	// Absolute Orientation does not work without lapack
}
#else // HAVE_LAPACK

void TestCovarianceAbsoluteOrientation()
{
	// do some iterations of random tests
	// testCovarianceAbsoluteOrientationRandom< float >( 10, 1e-1f );
	testCovarianceAbsoluteOrientationRandom< double >( 10, 1e-1 );
}

#endif // HAVE_LAPACK


#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utAlgorithm/ToolTip/Ransac.h>
#include <utAlgorithm/ToolTip/TipCalibration.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

template< typename T >
void testRobustTipCalibrationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -.5, .5 );
	typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.005, 0.005 ); // uniform noise for translation
	// typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.005 ); // gaussian noise for translation
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n = 10+(iRun%91);//( Random::distribute_uniform< std::size_t >( 4, 30 ) );
	
		//make up a random pose
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		const Pose pose( q, t );
		
		// choose some random position in space( the tip is orbiting)
		const Vector< T, 3 > origin = randVector();
		
		// generate some storage
		std::vector< Pose > poses;
		poses.reserve( n );
		
		std::vector< Pose > noisyPoses;
		noisyPoses.reserve( n );
		
		// generate some random poses via random rotations
		for( std::size_t i = 0; i<n; ++i )
		{
			const Pose newPose = Pose( randQuat(), origin ) * pose;
			poses.push_back( newPose );
			
			const T rotEps( 0.015 );
			const Pose noisyPose( Quaternion( newPose.rotation().x() + Random::distribute_uniform< T >( -rotEps, rotEps )
				, newPose.rotation().y() + Random::distribute_uniform< T >( -rotEps, rotEps )
				, newPose.rotation().z() + Random::distribute_uniform< T >( -rotEps, rotEps )
				, newPose.rotation().w() + Random::distribute_uniform< T >( -rotEps, rotEps ) ).normalize()
				, newPose.translation() + randPositionNoise() );
			
			noisyPoses.push_back( noisyPose );
		}
		
		// now produce some (10%) outlier
		const std::size_t outlier( n/10 );
		for( std::size_t i = 0; i<outlier; ++i )
		{
			const std::size_t index = Random::distribute_uniform< std::size_t >( 0, n-1 ) ;
			noisyPoses[ index ] = Pose( randQuat(), randVector() );
		}
		
		
		// set reasonable? ransac parameters
		const T threshold( 0.05 ); // <- this works quite well here, since there are some big outlier, might be differnet in other scenarios
		const std::size_t minSetSize( 3 ); // <- the algorihms needs at least 3 poses to estimate a result (minimum consensus set)
		const T percentOutlier( 0.2 ); // <- should be a reasonable number, mayber dervided empirically
		const T sucessProbability( 0.99 ); // <- usually you want to have this parameter quite high (~95 or higher)
		const Ubitrack::Math::Optimization::RansacParameter< T > params( threshold, minSetSize, n, percentOutlier, sucessProbability );
		//BOOST_MESSAGE( "Starting with at least " << params.nMaxIterations << " iterations for " << params.nMinInlier << " expected inlier of " << n << " poses."  );
		
		// do some robust estimation now
		Vector< T, 3 > pTool2Tip;
		Vector< T, 3 > pWolrd2Tip;
		if( !Ubitrack::Algorithm::ToolTip::estimatePosition3D_6D( pWolrd2Tip, noisyPoses, pTool2Tip, params ) )
		{
			// BOOST_MESSAGE( "robust tooltip calibration from " << n << " poses did not converge successfully with " << outlier << " outlier within " << params.nMaxIterations << " iterations." );
			continue;
		}
		
		// perform a common estimation for comparison
		Vector< T, 3 > pTool2Tip2;
		Vector< T, 3 > pWolrd2Tip2;
		if( !Ubitrack::Algorithm::ToolTip::estimatePosition3D_6D( pWolrd2Tip2, noisyPoses, pTool2Tip2 ) )
		{
			// BOOST_MESSAGE( "tooltip calibration did not converge successfully." );
			continue;
		}
		
		const std::pair< T, T > err1 = Ubitrack::Algorithm::ToolTip::estimatePosition3DError_6D( pWolrd2Tip, poses, pTool2Tip );
		const std::pair< T, T > err2 = Ubitrack::Algorithm::ToolTip::estimatePosition3DError_6D( pWolrd2Tip2, poses, pTool2Tip2 );
		
		BOOST_CHECK_MESSAGE( err1.first <= err2.first, "\nRobust tooltip calibration from " << n << " poses and " << outlier << " outlier resulted in a worse ERROR:\n" << err2.first << " (mean) +-" << err2.second << " (expected)\n" << err1.first << " (mean) +-" << err1.second << " (estimated)\n" );
		// BOOST_CHECK_MESSAGE( err1.first > err2.first, "\nRobust tooltip calibration using " << n << " poses resulted in SUCCESS:\n" << err2.first << ", " << err2.second << " (expected)\n" << err1.first << ", " << err1.second << " (estimated)\n" );
	}
}

void TestRobustTipCalibration()
{
	testRobustTipCalibrationRandom< double >( 1000, 1e-6 );
}


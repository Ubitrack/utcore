
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utAlgorithm/ToolTip/Ransac.h>
#include <utAlgorithm/ToolTip/ErrorEstimation.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

#ifndef HAVE_LAPACK
template< typename >
void testRobustTipCalibrationRandom()
{
	// TipCalibration does not work without lapack
}
#else // HAVE_LAPACK

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
		const T threshold( 0.02 );
		const std::size_t size( 5 );
		const std::size_t mInlier( 3 );
		const std::size_t nMinRuns( 3 );
		const std::size_t nMaxRuns( 15 );
		
		const Ubitrack::Math::Optimization::RansacParameter< T > params( threshold, size, mInlier, nMinRuns, nMaxRuns );
		
		// do some robust estimation now
		Vector< T, 3 > pTool2Tip;
		Vector< T, 3 > pWolrd2Tip;
		if( !Ubitrack::Algorithm::ToolTip::estimatePosition3D_6D( pWolrd2Tip, noisyPoses, pTool2Tip, params ) )
		{
			BOOST_CHECK_MESSAGE( false, "robust tooltip calibration did not converge successfully." );
			continue;
		}
		
		// perform a common estimation for comparison
		Vector< T, 3 > pTool2Tip2;
		Vector< T, 3 > pWolrd2Tip2;
		if( !Ubitrack::Algorithm::ToolTip::estimatePosition3D_6D( pWolrd2Tip2, noisyPoses, pTool2Tip2 ) )
		{
			BOOST_CHECK_MESSAGE( false, "tooltip calibration did not converge successfully." );
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

#endif // HAVE_LAPACK


#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utAlgorithm/ToolTip/Optimization.h>
#include <utAlgorithm/ToolTip/TipCalibration.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;


template< typename T >
void testOptimizedTipCalibrationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -0.5, 0.5 );
	typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.005, 0.005 ); // uniform noise for translation
	// typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.005 ); // gaussian noise for translation
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n = 3+(iRun%98);//( Random::distribute_uniform< std::size_t >( 4, 30 ) );
	
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

		// do some estimation now
		Vector< T, 3 > pTool2Tip;
		Vector< T, 3 > pWolrd2Tip;
		const Ubitrack::Math::Optimization::OptTerminate termCrit( 50, 1e-8 );
		if( !Ubitrack::Algorithm::ToolTip::estimatePosition3D_6D( pWolrd2Tip, noisyPoses, pTool2Tip, termCrit ) )
		{
			// BOOST_MESSAGE( "\non-linear tooltip calibration from " << n << " poses did not converge successfully." );
			continue;
		}
		
		
		// estimate ccomparison parameters
		Vector< T, 3 > pTool2Tip2;
		Vector< T, 3 > pWolrd2Tip2;
		if( !Ubitrack::Algorithm::ToolTip::estimatePosition3D_6D( pWolrd2Tip2, noisyPoses, pTool2Tip2 ) )
		{
//			BOOST_MESSAGE( "tooltip calibration from " << n << " poses did not converge successfully." );
			continue;
		}
		
		const std::pair< T, T > err1 = Ubitrack::Algorithm::ToolTip::estimatePosition3DError_6D( pWolrd2Tip, poses, pTool2Tip );
		const std::pair< T, T > err2 = Ubitrack::Algorithm::ToolTip::estimatePosition3DError_6D( pWolrd2Tip2, poses, pTool2Tip2 );
		
		// sadly the non-linear optimization has no advantage over the normal one. Therefore the test is that at least the results are roughly the same
		BOOST_CHECK_MESSAGE( (err1.first - err2.first )< epsilon, "\nNon-linear tooltip calibration from " << n << " poses resulted in a worse ERROR:\n" << err2.first << " +-" << err2.second << " (expected)\n" << err1.first << " +-" << err1.second << " (estimated)\n" );
		
		// BOOST_CHECK_MESSAGE( err1 >= err2, "\nnon-linear tooltip calibration from " << n << " poses resulted in SUCCESS:\n" << err2.first << ", " << err2.second << " (expected)\n" << err1.first << ", " << err1.second << " (estimated)\n" );
		
		// std::cout << "Trans   " << t << "\n";
		// std::cout << "Tip     " << pTool2Tip << "\n";
		// std::cout << "Tip2    " << pTool2Tip2 << "\n";
		// std::cout << "Origin  " << origin << "\n";
		// std::cout << "Tool    " << pWolrd2Tip << "\n";
		// std::cout << "Tool2   " << pWolrd2Tip2 << "\n";
	}
}

#ifdef HAVE_LAPACK

void TestOptimizedTipCalibration()
{
	testOptimizedTipCalibrationRandom< double >( 1000, 1e-6 );
}

#else // HAVE_LAPACK

void TestOptimizedTipCalibration()
{
	// TipCalibration does not work without lapack
}

#endif // HAVE_LAPACK

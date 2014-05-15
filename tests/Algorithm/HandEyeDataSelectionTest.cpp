
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>

#include <utAlgorithm/HandEye/DataSelection.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>

#include <algorithm>

#include "../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;


template< typename T >
void TestDataSelection( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -10., 10. );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n( Random::distribute_uniform< std::size_t >( 100, 200 ) );

		//generate a radnom pose
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		const Pose pose( q, t );
		
		//set first target frame
		std::vector< Pose > rightFrame;
		rightFrame.reserve( n );
		
		//set the second target frame
		std::vector< Pose > leftFrame;
		leftFrame.reserve( n );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			Quaternion q1 = randQuat();
			Vector< T, 3 > t1 = randVector();	
			Pose p1( q1, t1 );
			rightFrame.push_back( p1 );
			leftFrame.push_back( ~(pose * p1) );
		}

		// do some data Selection on poses
		const std::size_t diffs = (n*(n-1))/2;
		std::size_t selection( std::min< std::size_t >( diffs * 0.1, 1000u ) );
		std::vector< Pose > selectedPosesLeft;
		// selectedPosesLeft.reserve( selection );
		std::vector< Pose > selectedPosesRight;
		// selectedPosesRight.reserve( selection );
		
		Ubitrack::Algorithm::HandEye::select_6DPoses( leftFrame, rightFrame, selection, selectedPosesLeft, selectedPosesRight );
		// Ubitrack::Algorithm::select_6DPoses( leftFrame, rightFrame, selection, std::back_inserter( selectedPosesLeft ), std::back_inserter( selectedPosesRight ) );
		
		
		
		
		// if( b_done )
		// {
			// // calculate some errors
			// const T rotDiff = quaternionDiff( estimatedPose.rotation(), q );
			// const T posDiff = vectorDiff( estimatedPose.translation(), t );
			
			// // check if pose is better than before (only for valid results)
			// BOOST_CHECK_MESSAGE( rotDiff < epsilon, "\nEstimated rotation from " << n << " poses resulted in error " << rotDiff << " :\n" << q << " (expected)\n" << estimatedPose.rotation() << " (estimated)\n" );
			// BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nEstimated position from " << n << " poses resulted in error " << posDiff << " :\n" << t << " (expected)\n" << estimatedPose.translation() << " (estimated)\n" );
		// }
		// BOOST_WARN_MESSAGE( b_done, "Algorithm did not successfully estimate a result from " << n << " poses." );
	}
	BOOST_CHECK( 0 == 0 );
}

void TestHandEyeDataSelection()
{
	TestDataSelection< double >( 10, 1e-6 );
}

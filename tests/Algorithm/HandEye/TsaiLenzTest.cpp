
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/MatrixOperations.h>
#include <utAlgorithm/PoseEstimation6D6D/TsaiLenz.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

#ifndef HAVE_LAPACK
void TestHandEye()
{
	// HandyEye does not work without lapack
}
#else // HAVE_LAPACK

template< typename T >
void testHandEyeMatrixRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -100, 100 );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n( Random::distribute_uniform< std::size_t >( 4, 30 ) );

		//set first target frame
		std::vector< Matrix< T, 4, 4 > > rightFrame;
		rightFrame.reserve( n );
		for( std::size_t i = 0; i<n; ++i )
		{
			Quaternion q1 = randQuat();
			Vector< T, 3 > t1 = randVector();	
			Matrix< T, 4, 4 > mat1( q1, t1 );
			rightFrame.push_back( mat1 );
		}
		
		// produce other target frame
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		const Matrix< T, 4, 4 > mat( q, t );
		
		std::vector< Matrix< T, 4, 4 > > leftFrame;
		leftFrame.reserve( n );
		for( std::size_t i = 0; i<n; ++i )
		{
			// invert the Matrix
			const Matrix< T, 4, 4 > mTmp = boost::numeric::ublas::prod( mat, rightFrame[ i ] );
			leftFrame.push_back( invert_matrix( mTmp ) );
		}

		// do some estimation now
		const Pose estimatedPose = Ubitrack::Algorithm::PoseEstimation6D6D::performHandEyeCalibration ( leftFrame, rightFrame, true );
		
		// calculate some errors
		const T rotDiff = quaternionDiff( estimatedPose.rotation(), q );
		const T posDiff = vectorDiff( estimatedPose.translation(), t );
		// if( b_done )
		{
			// check if pose is better than before (only for valid results)
			BOOST_CHECK_MESSAGE( rotDiff < epsilon, "\nEstimated rotation from " << n << " poses resulted in error " << rotDiff << " :\n" << q << " (expected)\n" << estimatedPose.rotation() << " (estimated)\n" );
			BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nEstimated position from " << n << " poses resulted in error " << posDiff << " :\n" << t << " (expected)\n" << estimatedPose.translation() << " (estimated)\n" );
		}
		// BOOST_WARN_MESSAGE( b_done, "Algorithm did not successfully estimate a result with " << n 
			// << " points.\nRemaining difference in rotation " << rotDiff << ", difference in translation " << posDiff << "." );
	}
}

template< typename T >
void testHandEyePoseRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -10., 10. );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n( Random::distribute_uniform< std::size_t >( 4, 30 ) );

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

		// do some estimation now
		const Pose estimatedPose = Ubitrack::Algorithm::PoseEstimation6D6D::performHandEyeCalibration ( leftFrame, rightFrame, true );
		
		// calculate some errors
		const T rotDiff = quaternionDiff( estimatedPose.rotation(), q );
		const T posDiff = vectorDiff( estimatedPose.translation(), t );
		// if( b_done )
		{
			// check if pose is better than before (only for valid results)
			BOOST_CHECK_MESSAGE( rotDiff < epsilon, "\nEstimated rotation from " << n << " poses resulted in error " << rotDiff << " :\n" << q << " (expected)\n" << estimatedPose.rotation() << " (estimated)\n" );
			BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nEstimated position from " << n << " poses resulted in error " << posDiff << " :\n" << t << " (expected)\n" << estimatedPose.translation() << " (estimated)\n" );
		}
		// BOOST_WARN_MESSAGE( b_done, "Algorithm did not successfully estimate a result with " << n 
			// << " points.\nRemaining difference in rotation " << rotDiff << ", difference in translation " << posDiff << "." );
	}
}

void TestTsaiLenzHandEye()
{
	testHandEyeMatrixRandom< float >( 100, 1e-2f );
	testHandEyeMatrixRandom< double >( 100, 1e-6 );
	testHandEyePoseRandom< double >( 100, 1e-6 );
}

#endif // HAVE_LAPACK

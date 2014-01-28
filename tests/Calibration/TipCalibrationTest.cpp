
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utCalibration/TipCalibration.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

#ifndef HAVE_LAPACK
template< typename >
void testTipCalibrationRandom()
{
	// TipCalibration does not work without lapack
}
#else // HAVE_LAPACK

template< typename T >
void testTipCalibrationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -10., 10. );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n( Random::distribute_uniform< std::size_t >( 4, 30 ) );
		
		std::vector< Pose > poses;
		poses.reserve( n );

		//make up a random pose
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		const Pose pose( q, t );
		
		// generate some random poses via random rotations
		for( std::size_t i = 0; i<n; ++i )
		{
			const Pose newPose = Pose( randQuat(), Vector< T, 3 >( 0, 0 , 0 ) ) * pose;
			poses.push_back( newPose );
		}

		// do some estimation now
		Vector< T, 3 > pm;
		Vector< T, 3 > pw;
		Ubitrack::Calibration::tipCalibration( poses, pm, pw );
				
		for( std::size_t i = 0; i<n; ++i )
		{
			Vector< T, 3 > tip = poses[ i ] * pm;
			const T posDiff = vectorDistance( tip, pw );
			
			BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nEstimated tip position from " << n << " poses resulted in an error for one of the results" << posDiff << " :\n" << tip << " (expected)\n" << pw << " (estimated)\n" );
		}
	}
}

void TestTipCalibration()
{
	testTipCalibrationRandom< double >( 10000, 1e-6 );
}

#endif // HAVE_LAPACK

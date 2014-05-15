
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utAlgorithm/ToolTip/TipCalibration.h>
#include <utAlgorithm/ToolTip/ErrorEstimation.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

template< typename T >
void testTipCalibrationRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -.5, .5 );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n = 3+(iRun%98);
		
		std::vector< Pose > poses;
		poses.reserve( n );

		//make up a random pose
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		const Pose pose( q, t );
		
		// choose some random position in space( the tip is orbiting)
		const Vector< T, 3 > origin = randVector();
		
		// generate some random poses via random rotations
		for( std::size_t i = 0; i<n; ++i )
		{
			const Pose newPose = Pose( randQuat(), origin ) * pose;
			poses.push_back( newPose );
		}

		// do some estimation now
		Vector< T, 3 > pTool2Tip;
		Vector< T, 3 > pWolrd2Tip;
		if( !Ubitrack::Algorithm::ToolTip::estimatePosition3D_6D( pWolrd2Tip, poses, pTool2Tip ) )
		{
			BOOST_MESSAGE( "tooltip calibration from " << n << " poses did not converge successfully." );
			continue;
		}
		const std::pair< T, T > err = Ubitrack::Algorithm::ToolTip::estimatePosition3DError_6D( pWolrd2Tip, poses, pTool2Tip );
		
		BOOST_CHECK_MESSAGE( (err.first<epsilon), "\nTooltip calibration from " << n << " poses resulted in an estimated average ERROR:\n" << err.first << "(mean) +-" << err.second );
		
		// std::cout << "Origin  " << origin << "\n";
		// std::cout << "Tip     " << pTool2Tip << "\n";
		// std::cout << "Tool    " << pWolrd2Tip << "\n";
		
		// for( std::size_t i = 0; i<n; ++i )
		// {
			// const Vector< T, 3 > tip = poses[ i ] * pTool2Tip;
			// const T posDiff = vectorDistance( tip, pWolrd2Tip );
			
			// BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nEstimated tip position from " << n << " poses resulted in an error for one of the results " << posDiff << " :\n" << tip << " (expected)\n" << pWolrd2Tip << " (estimated)\n" );
		// }
	}
}

#ifdef HAVE_LAPACK

void TestTipCalibration()
{
	testTipCalibrationRandom< double >( 10000, 1e-6 );
}

#else // HAVE_LAPACK

template< typename >
void testTipCalibrationRandom()
{
	// TipCalibration does not work without lapack
}

#endif // HAVE_LAPACK

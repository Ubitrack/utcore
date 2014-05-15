
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointTransformation.h>
#include <utAlgorithm/PoseEstimation3D3D/AbsoluteOrientation.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

template< typename T >
void testRotation3DRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -100, 100 );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n( Random::distribute_uniform< std::size_t >( 3, 30 ) );

		std::vector< Vector< T, 3 > > rightFrame;
		rightFrame.reserve( n );
		std::generate_n ( std::back_inserter( rightFrame ), n,  randVector );
		
		
		Ubitrack::Math::Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		Matrix< T, 3, 4 > trafo( q, t );
		
		std::vector< Vector< T, 3 > > leftFrame;
		leftFrame.reserve( n );
		Geometry::transform_points( trafo, rightFrame.begin(), rightFrame.end(), std::back_inserter( leftFrame ) );

		// do some estimation now
		Ubitrack::Math::Quaternion estimatedQuat;
		const bool b_done = Ubitrack::Algorithm::PoseEstimation3D3D::estimateRotation_3D3D( leftFrame, estimatedQuat, rightFrame );
		
		// calculate some errors
		const T rotDiff = quaternionDiff( estimatedQuat, q );
		if( b_done )
		{
			// check if pose is better than before (only for valid results)
			BOOST_CHECK_MESSAGE( rotDiff < epsilon, "\nCompare rotation result (expected vs. estimated) using " << n << " points:\n" << q << " " << estimatedQuat );
		}
		BOOST_WARN_MESSAGE( b_done, "Algorithm did not successfully estimate a result with " << n 
			<< " points.\nRemaining difference in rotation " << rotDiff << "." );
	}
	
}

#ifndef HAVE_LAPACK

void TestAbsOrientRotation3D()
{
	// Absolute Orientation does not work without lapack
}

#else // HAVE_LAPACK

void TestAbsOrientRotation3D()
{
	// do some iterations of random tests
	testRotation3DRandom< float >( 10000, 1e-2f );
	testRotation3DRandom< double >( 10000, 1e-6 );
}

#endif // HAVE_LAPACK

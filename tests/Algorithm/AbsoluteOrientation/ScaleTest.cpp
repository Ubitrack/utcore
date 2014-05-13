
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointTransformation.h>
#include <utAlgorithm/AbsoluteOrientation/AbsoluteOrientation.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

#ifndef HAVE_LAPACK
void TestAbsOrientScale()
{
	// Absolute Orientation does not work without lapack
}
#else // HAVE_LAPACK


template< typename T >
void testScaleRandom( const std::size_t n_runs, const T epsilon )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -100, 100 );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n( Random::distribute_uniform< std::size_t >( 3, 30 ) );
		const T scale( Random::distribute_uniform< T >( 0.001, 100. ) );

		std::vector< Vector< T, 3 > > rightFrame;
		rightFrame.reserve( n );
		std::generate_n ( std::back_inserter( rightFrame ), n,  randVector );
		
		
		Ubitrack::Math::Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		Matrix< T, 3, 4 > poseMat( q, t );
		
		std::vector< Vector< T, 3 > > leftFrame;
		leftFrame.reserve( n );
		Geometry::transform_points( poseMat, rightFrame.begin(), rightFrame.end(), std::back_inserter( leftFrame ) );
		
		Matrix< T, 3, 3 > scaleMat = Matrix< T, 3, 3 >::identity();
		scaleMat( 0, 0 ) = scale;
		scaleMat( 1, 1 ) = scale;
		scaleMat( 2, 2 ) = scale;
		
		
		Geometry::transform_points( scaleMat, leftFrame.begin(), leftFrame.end(), leftFrame.begin() );

		// do some estimation now
		const T scale_e = Ubitrack::Algorithm::AbsoluteOrientation::estimateScale_3D3D( leftFrame, rightFrame );
		
		const T scale_diff( fabs(scale_e - scale) );
		// calculate some error
		BOOST_WARN_MESSAGE( scale_diff < epsilon, "Algorithm did not successfully estimate a result with " << n 
			<< " points.\nRemaining difference in scale estimation " << scale_diff << "." );
	}
	
}

void TestAbsOrientScale()
{
	// do some iterations of random tests
	testScaleRandom< float >( 10, 1e-2f );
	// testScaleRandom< double >( 10, 1e-6 );
}

#endif // HAVE_LAPACK

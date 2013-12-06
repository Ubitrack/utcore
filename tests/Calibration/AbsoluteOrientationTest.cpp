
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Functors/VectorFunctors.h>
#include <utCalibration/AbsoluteOrientation.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>

#include <utUtil/Exception.h>
#include "../tools.h"


#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace Ubitrack::Math;

#ifndef HAVE_LAPACK
void TestAbsoluteOrientation()
{
	// Absolute Orientation does not work without lapack
}
#else // HAVE_LAPACK

void fillDemoVectorsDeterministic( Vector< double, 3 >* left, Vector< double, 3 >* right, Quaternion q, Vector< double, 3 > t )
{
	left[0][0] = 1.0;
	left[0][1] = 0.0;
	left[0][2] = 0.0;

	left[1][0] = 0.0;
	left[1][1] = 1.0;
	left[1][2] = 0.0;

	left[2][0] = 0.0;
	left[2][1] = 0.0;
	left[2][2] = 1.0;

	left[3][0] = 0.5;
	left[3][1] = 0.5;
	left[3][2] = 0.7;

	right[0] = q*left[0]+t;
	right[1] = q*left[1]+t;
	right[2] = q*left[2]+t;
	right[3] = q*left[3]+t;
}

void testAbsoluteOrientationDeterministic()
{
	Vector< double, 3 > axis ( 1.0, 1.0, 1.5 );
	Quaternion q ( axis, M_PI/6.0 );
	Vector< double, 3 > t ( -1.0, 3.0, 2.5 );

	Vector< double, 3 > left[4];
	Vector< double, 3 > right[4];
	fillDemoVectorsDeterministic ( left, right, q, t );

	Pose p = Ubitrack::Calibration::calculateAbsoluteOrientation ( &left[0], &left[4], &right[0], &right[4] );

	BOOST_CHECK_SMALL ( vectorDiff ( p.translation(), t ), 10e-6 );
	BOOST_CHECK_SMALL ( quaternionDiff ( p.rotation(), q ), 10e-6 );
}

template< typename T >
void testAbsoluteOrientationRandom( const std::size_t n_runs, const T epsilon )
{

	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -100, 100 );
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		const std::size_t n( Random::distribute_uniform< std::size_t >( 4, 30 ) );

		std::vector< Vector< T, 3 > > leftFrame;
		leftFrame.reserve( n );
		std::generate_n ( std::back_inserter( leftFrame ), n,  randVector );
		
		
		Quaternion q = randQuat();
		Vector< T, 3 > t = randVector();
		
		std::vector< Vector< T, 3 > > rightFrame;
		rightFrame.reserve( n );
		std::transform( leftFrame.begin(), leftFrame.end(), std::back_inserter( rightFrame ), Functors::TransformVector< T >( q, t ) );

		Pose p = Ubitrack::Calibration::calculateAbsoluteOrientation ( leftFrame.begin(), leftFrame.end(), rightFrame.begin(), rightFrame.end() );
		BOOST_CHECK_SMALL ( vectorDiff ( p.translation(), t ), epsilon );
		BOOST_CHECK_SMALL ( quaternionDiff ( p.rotation(), q ), epsilon );
	}
	
}

void TestAbsoluteOrientation()
{
	// first do a deterministic test
	testAbsoluteOrientationDeterministic();
	
	// do some iterations of random tests
	testAbsoluteOrientationRandom< double >( 10000, 10e-6 );
}

#endif // HAVE_LAPACK

#define _USE_MATH_DEFINES

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <Math/Matrix.h>
#include <Math/Pose.h>
#include <Util/Exception.h>
#include <utCalibration/AbsoluteOrientation.h>
#include "../tools.h"
#include <cmath>

using namespace Ubitrack;

#ifndef HAVE_LAPACK
void TestAbsoluteOrientation()
{
	// Absolute Orientation does not work without lapack
}
#else // HAVE_LAPACK

void fillDemoVectorsDeterministic( Math::Vector<3>* left, Math::Vector<3>* right, Math::Quaternion q, Math::Vector< 3 > t )
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

void doDeterministicTest()
{
	Math::Vector<3> axis ( 1.0, 1.0, 1.5 );
	Math::Quaternion q ( axis, M_PI/6.0 );
	Math::Vector<3> t ( -1.0, 3.0, 2.5 );

	Math::Vector< 3 > left[4];
	Math::Vector< 3 > right[4];
	fillDemoVectorsDeterministic ( left, right, q, t );

	Math::Pose p = Calibration::calculateAbsoluteOrientation ( &left[0], &left[4], &right[0], &right[4] );

	BOOST_CHECK_SMALL ( vectorDiff ( p.translation(), t ), 10e-6 );
	BOOST_CHECK_SMALL ( quaternionDiff ( p.rotation(), q ), 10e-6 );
}

void doRandomizedTest()
{
	const int nPoints = 20;

	std::vector< Math::Vector< 3 > > leftFrame;
	std::vector< Math::Vector< 3 > > rightFrame;

	double magnification = random (0.1, 10.0);

	for (int j=0; j < nPoints; ++j)
	{
		leftFrame.push_back (randomVector<3, double>()*magnification);
	}

	Math::Quaternion q = randomQuaternion();
	Math::Vector < 3 > t = randomVector<3, double>();

	for (int j=0; j < nPoints; ++j)
	{
		rightFrame.push_back ( q*leftFrame[j]+t );
	}

	Math::Pose p = Calibration::calculateAbsoluteOrientation ( leftFrame.begin(), leftFrame.end(), rightFrame.begin(), rightFrame.end() );
	BOOST_CHECK_SMALL ( vectorDiff ( p.translation(), t ), 10e-6 );
	BOOST_CHECK_SMALL ( quaternionDiff ( p.rotation(), q ), 10e-6 );
}

void TestAbsoluteOrientation()
{
	// first do a deterministic test
	doDeterministicTest();

	// do some iterations of random tests
	const int nIterations = 100;
 	for (int i=0; i < nIterations; ++i)
 	{
		doRandomizedTest();
	}
}

#endif // HAVE_LAPACK

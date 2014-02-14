#define _USE_MATH_DEFINES

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <utUtil/Exception.h>
#include "../tools.h"

#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
//#include <utMath/Functors/VectorFunctors.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>

#include <utCalibration/Correlation.h>


using namespace Ubitrack::Math;

#ifndef HAVE_LAPACK
void TestCorrelation()
{
}
#else // HAVE_LAPACK

std::vector<double>* generateEmptyPositionSequence() {
	std::vector<double> *v = new std::vector<double>();
	return v;
}

double randDouble() {
	return Random::distribute_uniform<double>(0.0f, 1.0f);
}

std::vector<double>* generateRandomPositionSequence() {
	std::vector<double> *v = new std::vector<double>();
	v->reserve(100);

	std::generate_n( std::back_inserter(*v), 100, randDouble );

	return v;
}


void TestCorrelation()
{
	{
		std::vector<double> *v1 = generateEmptyPositionSequence();
		BOOST_CHECK( Ubitrack::Calibration::computeCorrelation(*v1, *v1) == 1.0f );
		delete v1;
	}

	{
		std::vector<double> *v1 = generateRandomPositionSequence();
		BOOST_CHECK( Ubitrack::Calibration::computeCorrelation(*v1, *v1) == 1.0f );
		delete v1;
	}

	{
		std::vector<double> *v1 = generateRandomPositionSequence();
		std::vector<double> *v2 = generateRandomPositionSequence();

		BOOST_CHECK( Ubitrack::Calibration::computeCorrelation(*v1, *v2) != 1.0f );
		BOOST_CHECK( Ubitrack::Calibration::computeCorrelation(*v1, *v2) < 1.0f );
		delete v1;
		delete v2;
	}
}

#endif // HAVE_LAPACK

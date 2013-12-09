#define _USE_MATH_DEFINES

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <utUtil/Exception.h>
#include "../tools.h"

#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Functors/VectorFunctors.h>

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

std::vector<Vector3d>* generateEmptyPositionSequence() {
	std::vector<Vector3d> *v = new std::vector<Vector3d>();
	return v;
}

std::vector<Vector3d>* generateRandomPositionSequence() {
	std::vector<Vector3d> *v = new std::vector<Vector3d>();

	for (int i=0; i<100; i++) {
		
	}

	return v;
}


void TestCorrelation()
{
	std::vector<Vector3d> *v1 = generateEmptyPositionSequence();

	BOOST_CHECK( Ubitrack::Calibration::computeCorrelation(*v1, *v1) == 1.0f );

}

#endif // HAVE_LAPACK


#include <boost/test/unit_test.hpp>
#include <utUtil/Logging.h>
#include "Math/MathTest.h"
#include "MathUtil/MathUtilTest.h"
#include "Geometry/GeometryTest.h"
#include "Stochastic/StochasticTest.h"
#include "Algorithm/AlgorithmTest.h"

using boost::unit_test::test_suite;

// this function replaces the C++ main function, which is implemented by BOOST
// run with --log_level=all to get some more output
test_suite* init_unit_test_suite( int, char* [] )
{
#ifndef NDEBUG
	// Ubitrack::Util::initLogging();
#endif
	
	// create a test suite
	test_suite* allTests = BOOST_TEST_SUITE( "utcore" );

	// this example will pass cause we know ahead of time number of expected failures
	allTests->add( new MathTest );
	allTests->add( new MathUtilTest );
	allTests->add( new GeometryTest );
	allTests->add( new StochasticTest );
	allTests->add( new AlgorithmTest );
	
	return allTests;
}


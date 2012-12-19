#include <boost/test/unit_test.hpp>
#include <utUtil/Logging.h>
#include "Math/MathTest.h"
#include "Calibration/CalibTest.h"

using boost::unit_test::test_suite;

// this function replaces the C++ main function, which is implemented by BOOST
// run with --log_level=all to get some more output
test_suite* init_unit_test_suite( int, char* [] )
{
	Ubitrack::Util::initLogging();
	
	// create a test suite
	test_suite* allTests = BOOST_TEST_SUITE( "Module utcore" );

	// this example will pass cause we know ahead of time number of expected failures
	allTests->add( new MathTest );
	allTests->add( new CalibrationTest );

	return allTests;
}


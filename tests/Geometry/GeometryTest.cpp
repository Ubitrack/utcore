#include "GeometryTest.h"

// declare external tests here, to save us some trivial header files
void TestPoints();
void TestConic();


GeometryTest::GeometryTest()
	: boost::unit_test::test_suite( "GeometryTests" )
{
	add( BOOST_TEST_CASE( &TestPoints ) );
	add( BOOST_TEST_CASE( &TestConic ) );
}

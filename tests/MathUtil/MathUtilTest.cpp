#include "MathUtilTest.h"

// declare external tests here, to save us some trivial header files
void TestRotationCast();

MathUtilTest::MathUtilTest()
	: boost::unit_test::test_suite( "MathUtilTests" )
{
	add( BOOST_TEST_CASE( &TestRotationCast ) );
}

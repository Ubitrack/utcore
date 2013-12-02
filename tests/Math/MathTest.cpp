#include "MathTest.h"

// declare external tests here, to save us some trivial header files
void TestConic();
void TestLapack();
void TestQuaternionConversion();


MathTest::MathTest()
	: boost::unit_test::test_suite( "Math test suite" )
{
	add( BOOST_TEST_CASE( &TestConic ) );
	add( BOOST_TEST_CASE( &TestLapack ) );
	add( BOOST_TEST_CASE( &TestQuaternionConversion ) );
}

#include "MathTest.h"

// declare external tests here, to save us some trivial header files
void TestBlas1();
void TestBlas2();
void TestVectorFunctions();
void TestLapack();


MathTest::MathTest()
	: boost::unit_test::test_suite( "MathTests" )
{
	add( BOOST_TEST_CASE( &TestBlas1 ) );
	add( BOOST_TEST_CASE( &TestBlas2 ) );
	add( BOOST_TEST_CASE( &TestVectorFunctions ) );
	add( BOOST_TEST_CASE( &TestLapack ) );
}

#include "StochasticTest.h"

// declare external tests here, to save us some trivial header files
void TestKMeans();



StochasticTest::StochasticTest()
	: boost::unit_test::test_suite( "Stochastic test suite" )
{
	add( BOOST_TEST_CASE( &TestKMeans ) );
}

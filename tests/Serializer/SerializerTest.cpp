#include "SerializerTest.h"

// declare external tests here, to save us some trivial header files
void TestBoostArchive();


SerializerTest::SerializerTest()
	: boost::unit_test::test_suite( "SerializerTests" )
{
	add( BOOST_TEST_CASE( &TestBoostArchive ) );
}

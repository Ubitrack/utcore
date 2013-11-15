#include "CalibTest.h"

// declare external tests here, to save us some trivial header files
void Test2D3DPoseEstimation();
void Test3DPointReconstruction();
void TestAbsoluteOrientation();
void TestBundleAdjustment();
void TestDecomposeProjection();
void TestFundamentalMatrix();
void TestHomographyDLT();
void TestProjectionDLT();
void TestSquareHomography();

CalibrationTest::CalibrationTest()
	: boost::unit_test::test_suite( "Calibration test suite" )
{
	add( BOOST_TEST_CASE( &Test2D3DPoseEstimation ) );
	add( BOOST_TEST_CASE( &Test3DPointReconstruction ) );
	add( BOOST_TEST_CASE( &TestAbsoluteOrientation ) );
	add( BOOST_TEST_CASE( &TestBundleAdjustment ) );
	add( BOOST_TEST_CASE( &TestDecomposeProjection ) );
	add( BOOST_TEST_CASE( &TestFundamentalMatrix ) );
	add( BOOST_TEST_CASE( &TestHomographyDLT ) );
	add( BOOST_TEST_CASE( &TestProjectionDLT ) );
	add( BOOST_TEST_CASE( &TestSquareHomography ) );
}


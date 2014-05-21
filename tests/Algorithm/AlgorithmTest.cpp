#include "AlgorithmTest.h"

// declare external tests here, to save us some trivial header files
// new tests
void TestTipCalibration();
void TestRobustTipCalibration();
void TestOptimizedTipCalibration();
void TestAbsOrientScale();
void TestAbsOrientRotation3D();
void TestAbsoluteOrientation();
void TestRobustAbsoluteOrientation();
void TestOptimizedAbsoluteOrientation();
void TestCovarianceAbsoluteOrientation();

// old tests..
void Test2D3DPoseEstimation();
void Test3DPointReconstruction();
void TestBundleAdjustment();
void TestDecomposeProjection();
void TestFundamentalMatrix();
void TestHomography();
void TestProjectionDLT();
void TestCorrelation();
void TestHandEye();
void TestDualHandEye();
void TestHandEyeDataSelection();

AlgorithmTest::AlgorithmTest()
	: boost::unit_test::test_suite( "AlgorithmTests" )
{
	add( BOOST_TEST_CASE( &TestTipCalibration ) );
	add( BOOST_TEST_CASE( &TestRobustTipCalibration ) );
	add( BOOST_TEST_CASE( &TestOptimizedTipCalibration ) );
	add( BOOST_TEST_CASE( &TestAbsOrientScale ) );
	add( BOOST_TEST_CASE( &TestAbsOrientRotation3D ) );
	add( BOOST_TEST_CASE( &TestAbsoluteOrientation ) );
	add( BOOST_TEST_CASE( &TestRobustAbsoluteOrientation ) );
	add( BOOST_TEST_CASE( &TestOptimizedAbsoluteOrientation ) );
	add( BOOST_TEST_CASE( &TestCovarianceAbsoluteOrientation ) );
	
	// old tests...
	add( BOOST_TEST_CASE( &Test2D3DPoseEstimation ) );
	add( BOOST_TEST_CASE( &Test3DPointReconstruction ) );
	add( BOOST_TEST_CASE( &TestBundleAdjustment ) );
	add( BOOST_TEST_CASE( &TestDecomposeProjection ) );
	add( BOOST_TEST_CASE( &TestFundamentalMatrix ) );
	add( BOOST_TEST_CASE( &TestHomography ) );
	add( BOOST_TEST_CASE( &TestProjectionDLT ) );
	add( BOOST_TEST_CASE( &TestHandEye ) );
	add( BOOST_TEST_CASE( &TestDualHandEye ) );
	add( BOOST_TEST_CASE( &TestHandEyeDataSelection ) );
	add( BOOST_TEST_CASE( &TestCorrelation ) );
	

}


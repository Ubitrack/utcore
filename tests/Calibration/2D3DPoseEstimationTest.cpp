#include <boost/test/unit_test.hpp>
#include <math.h>
#include <boost/test/floating_point_comparison.hpp>
#include <utCalibration/2D3DPoseEstimation.h>
#include "../tools.h"
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <iostream>

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;

void TestOptimizePose()
{
	for ( int iRun = 0; iRun < 10; iRun++ )
	{
		// random intrinsics matrix
		Matrix< 3, 3, double > cam( ublas::identity_matrix< double >( 3 ) );
		cam( 0, 0 ) = random( 200, 800 );
		cam( 1, 1 ) = random( 200, 800 );
		
		// random pose
		Quaternion rot( randomQuaternion() );
		Vector< 3, double > trans( randomVector< 3, double >() );
		trans( 2 ) = random( 1.0, 2.0 );
		
		// some random 3d points
		std::vector< Vector< 3, double > > p3D( 15 );
		for ( unsigned i = 0; i < p3D.size(); i++ )
			p3D[ i ] = randomVector< 3 >( 0.4 );
			
		// project to 2D points 
		std::vector< Vector< 2, double > > p2D( p3D.size() );
		for ( unsigned i = 0; i < p3D.size(); i++ )
		{
			Vector< 3, double > tmp( rot * p3D[ i ] + trans );
			tmp = ublas::prod( cam, tmp );
			p2D[ i ] = ublas::subrange( tmp, 0, 2 ) / tmp( 2 );
		}
		
		// add some noise to the pose
		Pose testPose( 
			Quaternion( rot.x() + random( -0.1, 0.1 ), rot.y() + random( -0.1, 0.1 ), 
				rot.z() + random( -0.1, 0.1 ), rot.w() + random( -0.1, 0.1 ) ),
			trans + randomVector< 3 >( 0.2 ) );
			
		Pose optimized( testPose );
			
		// run optimizePose
		Ubitrack::Calibration::optimizePose( optimized, p2D, p3D, cam );
		
		// check if pose is better than before
		BOOST_CHECK( quaternionDiff( testPose.rotation(), rot ) >= quaternionDiff( optimized.rotation(), rot ) );
		BOOST_CHECK( ublas::norm_2( testPose.translation() - trans ) >= ublas::norm_2( optimized.translation() - trans ) );
		BOOST_CHECK_SMALL( quaternionDiff( optimized.rotation(), rot ), 1e-3 );
		BOOST_CHECK_SMALL( ublas::norm_2( optimized.translation() - trans ), 1e-3 );
	}
}

void Test2D3DPoseEstimation()
{
	TestOptimizePose();
}


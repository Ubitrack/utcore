#include <boost/test/unit_test.hpp>
#include <math.h>
#include <boost/test/floating_point_comparison.hpp>
#include <utCalibration/2D3DPoseEstimation.h>
#include "../tools.h"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <utMath/RandomNumbers.h>
#include <utMath/Functors/UniformDistribution.h>

#include <iostream>
#include <algorithm>

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;

template< typename T >
void TestOptimizePose( const std::size_t n_runs, const T epsilon )
{

	Functors::uniform_quaternion<> randQuat;
	Functors::uniform_distribution< 3, T > randVector( -0.4, 0.4 );
	Functors::uniform_distribution< 3, T > randTranslation( -100, 100 );
	
	for ( int iRun = 0; iRun < n_runs; iRun++ )
	{
		// random intrinsics matrix
		Matrix< 3, 3, T > cam( ublas::identity_matrix< T >( 3 ) );
		cam( 0, 0 ) = distribute_uniform< T >( 200, 800 );
		cam( 1, 1 ) = distribute_uniform< T >( 200, 800 );
		
		// random pose
		Quaternion rot( randQuat( ) );
		Vector< 3, T > trans ( randTranslation() );
		trans( 2 ) = distribute_uniform< T >( 1, 10 );
		
		// some random 3d points
		const std::size_t n( distribute_uniform< std::size_t >( 5, 30 ) );
		std::vector< Ubitrack::Math::Vector< 3, T > > p3D;
		p3D.reserve( n );
		std::generate_n ( std::back_inserter( p3D ), n,  randVector );
		// std::copy( p3D.begin(), p3D.end(), std::ostream_iterator< Ubitrack::Math::Vector< 3 > > ( std::cout, ", ") );
			
		// project to 2D points 
		std::vector< Vector< 2, T > > p2D;
		p2D.reserve( n );
		for ( std::size_t i = 0; i < n; i++ )
		{
			Vector< 3, T > tmp( rot * p3D[ i ] + trans );
			tmp = ublas::prod( cam, tmp );
			p2D.push_back( ublas::subrange( tmp, 0, 2 ) / tmp( 2 ) );
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
		BOOST_CHECK_SMALL( quaternionDiff( optimized.rotation(), rot ), epsilon );
		BOOST_CHECK_SMALL( ublas::norm_2( optimized.translation() - trans ), epsilon );
	}
}

void Test2D3DPoseEstimation()
{
	TestOptimizePose< double >( 10, 1e-3 );
}


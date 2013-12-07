#include <boost/test/unit_test.hpp>
#include <math.h>
#include <boost/test/floating_point_comparison.hpp>
#include <utCalibration/2D3DPoseEstimation.h>
#include "../tools.h"
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Functors/VectorFunctors.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>



#include <iostream>
#include <algorithm>

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;

template< typename T >
void TestOptimizePose( const std::size_t n_runs, const T epsilon )
{

	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -0.5, 0.5 ); // 3d Points
	typename Random::Vector< T, 3 >::Uniform randTranslation( -100, 100 ); //translation
	typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.2 ); // translation gaussian noise
	// typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.3, 0.3 ); // translation uniform noise
	
	
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		// random intrinsics matrix
		Matrix< T, 3, 3 > cam( Matrix< T, 3, 3 >::identity() );
		cam( 0, 0 ) = Random::distribute_uniform< T >( 200, 800 );
		cam( 1, 1 ) = Random::distribute_uniform< T >( 200, 800 );
		
		// random pose
		Quaternion rot( randQuat( ) );
		Vector< T, 3 > trans ( randTranslation() );
		trans( 2 ) = Random::distribute_uniform< T >( 10, 100 );
		
		// some random 3d points
		// const std::size_t n( 15 );
		const std::size_t n( Random::distribute_uniform< std::size_t >( 10, 50 ) );
		
		
		std::vector< Ubitrack::Math::Vector< T, 3 > > p3D;
		p3D.reserve( n );
		std::generate_n ( std::back_inserter( p3D ), n,  randVector );
		// std::copy( p3D.begin(), p3D.end(), std::ostream_iterator< Ubitrack::Math::Vector< double, 3 > > ( std::cout, ", ") );
			
		// project to 2D points 
		Matrix< T, 3, 4 > proj( rot, trans );
		proj = boost::numeric::ublas::prod( cam, proj );
		std::vector< Vector< T, 2 > > p2D;
		p2D.reserve( n );
		std::transform( p3D.begin(), p3D.end(), std::back_inserter( p2D ), Functors::ProjectVector< T >( proj ) );
		// std::copy( p2D.begin(), p2D.end(), std::ostream_iterator< Ubitrack::Math::Vector< double, 2 > > ( std::cout, ", ") );
		
		// add some noise to the pose
		Pose testPose( 
			Quaternion( rot.x() + Random::distribute_uniform< T >( -0.1, 0.1 )
						, rot.y() + Random::distribute_uniform< T >( -0.1, 0.1 )
						, rot.z() + Random::distribute_uniform< T >( -0.1, 0.1 )
						, rot.w() + Random::distribute_uniform< T >( -0.1, 0.1 ) )
						, trans + randPositionNoise() );
			
		Pose optimized( testPose );

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
	TestOptimizePose< double >( 1000, 1e-3 );
	// TestOptimizePose< float >( 1000, 1e-2f );
}


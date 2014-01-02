

#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointProjection.h>
#include <utCalibration/2D3DPoseEstimation.h>
#include <utCalibration/2D3DPoseEstimationHager.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../tools.h"

#include <math.h>
#include <iostream>
#include <algorithm>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

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
		
		//generate the projection
		Matrix< T, 3, 4 > proj( rot, trans );
		proj = boost::numeric::ublas::prod( cam, proj );
		
		// some random 3d points
		// const std::size_t n( 15 );
		const std::size_t n( Random::distribute_uniform< std::size_t >( 10, 50 ) );
		
		
		std::vector< Ubitrack::Math::Vector< T, 3 > > p3D;
		p3D.reserve( n );
		std::generate_n ( std::back_inserter( p3D ), n,  randVector );
		// std::copy( p3D.begin(), p3D.end(), std::ostream_iterator< Ubitrack::Math::Vector< double, 3 > > ( std::cout, ", ") );
			
		// project to 2D points 		
		std::vector< Vector< T, 2 > > p2D;
		p2D.reserve( n );
		Geometry::project_points( proj, p3D.begin(), p3D.end(), std::back_inserter( p2D ) );
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
		const T rotDiff = quaternionDiff( optimized.rotation(), rot );
		const T posDiff = ublas::norm_2( optimized.translation() - trans );
		BOOST_WARN_SMALL( rotDiff, epsilon );
		BOOST_WARN_SMALL( posDiff, epsilon );
		BOOST_CHECK( quaternionDiff( testPose.rotation(), rot ) >= rotDiff );
		BOOST_CHECK( ublas::norm_2( testPose.translation() - trans ) >= posDiff );
	}
}

template< typename T >
void Test2D3DPoseEstimationGeneral( const std::size_t n_runs, const T epsilon )
{
	//tried to reassemble a typical use case
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randVector( -0.5, 0.5 ); // 3d Points
	typename Random::Vector< T, 3 >::Uniform randTranslation( -2, 2 ); //translation	
	
	std::size_t iter_count = 0;
	for ( std::size_t iRun = 0; iRun < n_runs; iRun++ )
	{
		// random pose
		Quaternion rot( randQuat( ) );
		Vector< T, 3 > trans ( randTranslation() );
		trans( 2 ) = Random::distribute_uniform< T >( 2, 5 );
		
		// projection
		Matrix< T, 3, 4 > proj( rot, trans );
		
		// some random 3d points
		// const std::size_t n( 5 );
		// attention:: small number of points leads to bigger errors
		// nevertheless it is still quite reliable
		const std::size_t n( Random::distribute_uniform< std::size_t >( 7, 50 ) );
		
		std::vector< Ubitrack::Math::Vector< T, 3 > > p3D;
		p3D.reserve( n );
		std::generate_n ( std::back_inserter( p3D ), n,  randVector );
		
		// project to 2D points 
		std::vector< Vector< T, 2 > > p2D;
		p2D.reserve( n );
		Geometry::project_points( proj, p3D.begin(), p3D.end(), std::back_inserter( p2D ) );

		//estimate the "unknown" pose
		Pose estimatedPose;
		T max_error( 1e-06 );
		std::size_t max_iterations( 100 );
		const bool b_done = Ubitrack::Calibration::estimatePose6D_2D3D( p2D, estimatedPose, p3D, max_iterations, max_error );
		iter_count += max_iterations;
		
		// estimate the differences
		const T rotDiff = quaternionDiff( estimatedPose.rotation(), rot );
		const T posDiff = ublas::norm_2( estimatedPose.translation() - trans );
		if( b_done )
		{
			// check if pose is better than before (only for valid results)
			BOOST_CHECK_MESSAGE( rotDiff < epsilon, "\nCompare rotation    result (expected vs. estimated) after " << max_iterations << " iterations using " << n << " points:\n" << Pose( proj ).rotation() << " " << estimatedPose.rotation() );
			BOOST_CHECK_MESSAGE( posDiff < epsilon, "\nCompare translation result (expected vs. estimated) after " << max_iterations << " iterations using " << n << " points:\n" << Pose( proj ).translation() << " " << estimatedPose.translation() );
		}
		BOOST_WARN_MESSAGE( b_done, "Algorithm did not converge after " << max_iterations << " iterations with " << n 
			<< " points.\nRemaining difference in rotation " << rotDiff << ", difference in translation " << posDiff << "." );
	}
	BOOST_MESSAGE( "Average number of iterations after " << n_runs << " runs: " << iter_count / n_runs );
}

void Test2D3DPoseEstimation()
{
	TestOptimizePose< float >( 1000, 1e-1f );
	TestOptimizePose< double >( 1000, 1e-3 );
	Test2D3DPoseEstimationGeneral< float >( 1000, 1e-01 );
	Test2D3DPoseEstimationGeneral< double >( 1000, 1e-01 );
}


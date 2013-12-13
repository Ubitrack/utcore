


#include <utMath/Pose.h>
#include <utMath/Scalar.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Functors/VectorFunctors.h> // e.g. Norm_2



#include <utMath/Random/Pose.h>
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include <utCalibration/BundleAdjustment.h> //no test yet, same as MultiCamerPoseOtpimization :(
//#include <utCalibration/MultipleCameraPoseOptimization.h>
#include "../tools.h"

#include <iostream>
#include <algorithm>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>



using namespace Ubitrack::Math;

namespace { //anonymous namespace

template< typename T >
struct isPixelNotWithinScreen
{
protected:
	const Ubitrack::Math::Vector< T, 2 > m_maxDimension;
	
public:
	isPixelNotWithinScreen( const Ubitrack::Math::Vector< T, 2 >& screenResolution )
		: m_maxDimension( screenResolution )
		{}
	
	bool operator()( Ubitrack::Math::Vector< T, 2 > &pixel ) const
	{
		// just check if pixel is within given screen resolution.
		if ( pixel[ 0 ] < 0 || pixel [ 1 ] < 0 || pixel[ 0 ] >= m_maxDimension[ 0 ] || pixel [ 1 ] >= m_maxDimension[ 1 ] )
			return true;
		return false;
	}
};

} // anonymous namesapce

template< typename T >
void TestMarkerBundleAdjustment( const std::size_t n_runs, const T epsilon )
{
	Vector< T, 2 > screenResolution( 640, 480 );
	
	
	// run bundleAdjustment test for several times
	for( std::size_t i( 0 ); i< n_runs; ++i )
	{
		//change here for a really big bundle adjustment:
		const std::size_t n_p3d = 10;//Random::distribute_uniform( 10, 15 ); 
		const std::size_t n_cams = Random::distribute_uniform( 3, 10 );
		
		// random intrinsics matrix, never changes, assume always the same camera
		Matrix< T, 3, 3 > cam( Matrix< T, 3, 3 >::identity() );
		cam( 0, 0 ) = Random::distribute_uniform< T >( 500, 800 );
		cam( 1, 1 ) = Random::distribute_uniform< T >( 500, 800 );
		cam( 0, 2 ) = -screenResolution[ 0 ] * 0.5;
		cam( 1, 2 ) = -screenResolution[ 1 ] * 0.5;
		cam( 2, 2 ) = -1;
			
		typename Random::Vector< T, 2 >::Normal randPixelNoise( 0, 0.25 ); // gaussian noise for 2d Pixels -> usually very below 1 pixel
		typename Random::Vector< T, 3 >::Uniform randVector( -10.0, 10.0 ); // 3d points, assume meters as unit
		typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.10 ); // gaussian noise for 3d points -> can be quite high
		// typename Random::Vector< T, 3 >::Uniform randPositionNoise( -0.25, 0.25 ); // uniform noise	for 3d points
		
		// random poses for the estimation
		typename Random::Quaternion< T >::Uniform randQuat;
		typename Random::Vector< T, 3 >::Uniform randTranslation( -10, 10 ); //translation
		//alternatively one could use a random pose:
		// Random::Pose< T >::Uniform( -100, 100 );
		
	
		// generate 3d points (ground truth) + add some noise to points
		std::vector< Vector< T, 3 > > points_3D;
		points_3D.reserve( n_p3d );
		std::generate_n ( std::back_inserter( points_3D ), n_p3d,  randVector );
		
		//generate some noise for the 3d points and add it
		std::vector< Vector< T, 3 > > noisy_points_3D;
		noisy_points_3D.reserve( n_p3d );
		std::generate_n ( std::back_inserter( noisy_points_3D ), n_p3d,  randPositionNoise );
		std::transform( points_3D.begin(), points_3D.end(), noisy_points_3D.begin(), noisy_points_3D.begin(), std::plus< Vector< T, 3 > >() );
		
		// generate random poses and project on image plane
		std::vector< Pose > extrinsics;
		std::vector< Matrix< T, 3, 3 > > intrinsics;
		std::vector< std::vector< Vector< T, 2 > > > observed_points_2D;
		
		do
		{
			Quaternion rot( randQuat( ) );
			Vector< T, 3 > trans ( randTranslation() );
			
			Pose pose6D( rot, trans );
			Matrix< T, 3, 4 > proj( rot, trans );
			proj = boost::numeric::ublas::prod( cam, proj );
			
			std::vector< Vector< T, 2 > > points_2D;
			points_2D.reserve( n_p3d );
			std::transform( points_3D.begin(), points_3D.end(), std::back_inserter( points_2D ), Functors::ProjectVector< T >( proj ) );
			
			//generate some noise for the 2d Points
			std::vector< Vector< T, 2 > > noisy_points_2D;
			noisy_points_2D.reserve( n_p3d );
			std::generate_n ( std::back_inserter( noisy_points_2D ), n_p3d,  randPixelNoise );
			std::transform( points_2D.begin(), points_2D.end(), noisy_points_2D.begin(), noisy_points_2D.begin(), std::plus< Vector< T, 2 > >() );
			
			// remove pixels outside the screen
			//noisy_points_2D.erase( std::remove_if( noisy_points_2D.begin(), noisy_points_2D.end(), isPixelNotWithinScreen< T >( screenResolution ) ), noisy_points_2D.end() ); 
			const std::size_t n_2d( noisy_points_2D.size() );
			
			if( n_2d > 4 ) //enough points are visible?
			{
				// Pose noisyPose( Quaternion( rot.x() + Random::distribute_uniform< T >( -0.1, 0.1 )
						// , rot.y() + Random::distribute_uniform< T >( -0.1, 0.1 )
						// , rot.z() + Random::distribute_uniform< T >( -0.1, 0.1 )
						// , rot.w() + Random::distribute_uniform< T >( -0.1, 0.1 ) )
						// , trans + randPositionNoise() );
						
				extrinsics.push_back( pose6D );
				intrinsics.push_back( cam );
				observed_points_2D.push_back( noisy_points_2D );
			}
		}
		while( extrinsics.size() < n_cams );
	
		
		
		std::vector< Ubitrack::Math::Vector< T, 3 > > vecBefore;
		vecBefore.reserve( n_p3d );
		std::transform( noisy_points_3D.begin(), noisy_points_3D.end(), points_3D.begin(), std::back_inserter( vecBefore ), std::minus< Ubitrack::Math::Vector< T, 3 > >() );
		
		std::vector< T > distsBefore;
		distsBefore.reserve( n_p3d );
		std::transform( vecBefore.begin(), vecBefore.end(), std::back_inserter( distsBefore ), Ubitrack::Math::Functors::Norm_2< T, 3 >() );
		const T distBefore = std::accumulate( distsBefore.begin(), distsBefore.end(), static_cast< T > ( 0 ) ) / n_p3d;
		
		
		Ubitrack::Calibration::simpleBundleAdjustment( observed_points_2D, intrinsics, extrinsics, noisy_points_3D );
		// call to templated function (does not link on windows):
		// Ubitrack::Calibration::simpleBundleAdjustment( observed_points_2D.begin(), observed_points_2D.end(), intrinsics.begin(), extrinsics.begin(), noisy_points_3D.begin(), noisy_points_3D.end() );
		
		
		
		std::vector< Ubitrack::Math::Vector< T, 3 > > vecAfter;
		vecAfter.reserve( n_p3d );
		std::transform( noisy_points_3D.begin(), noisy_points_3D.end(), points_3D.begin(), std::back_inserter( vecAfter ), std::minus< Ubitrack::Math::Vector< T, 3 > >() );
		std::vector< T > distsAfter;
		distsAfter.reserve( n_p3d );
		std::transform( vecAfter.begin(), vecAfter.end(), std::back_inserter( distsAfter ), Ubitrack::Math::Functors::Norm_2< T, 3 >() );
		const T distAfter = std::accumulate( distsAfter.begin(), distsAfter.end(), static_cast< T > ( 0 ) ) / n_p3d;
		
		
		BOOST_CHECK( distBefore >= distAfter );

		// // print out 3D error
		// std::cout << std::endl << "Distance before vs. after optimization: " << distBefore << " vs. " << distAfter << std::endl << std::endl;
		
		// // print out 3d values
		// std::cout << std::endl << std::endl;
		// std::copy( noisy_points_3D.begin(), noisy_points_3D.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 3 > > ( std::cout, ", ") );	
		
	
		
		// std::cout << std::endl << std::endl;
		// std::copy( extrinsics.begin(), extrinsics.end(), std::ostream_iterator< Ubitrack::Math::Pose > ( std::cout, ", ") );	
		// std::cout << std::endl;
		// std::copy( errPoses.begin(), errPoses.end(), std::ostream_iterator< Ubitrack::Math::ErrorPose > ( std::cout, ", ") );
		
		// check if pose is better than before
		// BOOST_CHECK( quaternionDiff( testPose.rotation(), rot ) >= quaternionDiff( optimized.rotation(), rot ) );
		// BOOST_CHECK( ublas::norm_2( testPose.translation() - trans ) >= ublas::norm_2( optimized.translation() - trans ) );
		// BOOST_CHECK_SMALL( quaternionDiff( optimized.rotation(), rot ), epsilon );
		// BOOST_CHECK_SMALL( ublas::norm_2( optimized.translation() - trans ), epsilon );
	}	
};

void TestBundleAdjustment()
{
	// attention works also with float now :)
	TestMarkerBundleAdjustment< double >( 10, 1e-3 );
	TestMarkerBundleAdjustment< float >( 10, 1e-3 );
}


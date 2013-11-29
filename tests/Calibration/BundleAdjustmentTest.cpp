


#include <utMath/Pose.h>
#include <utMath/Scalar.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Functors/VectorFunctors.h>

#include <utMath/Random/Pose.h>
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
//#include <utCalibration/BundleAdjustment.h> no test yet
#include <utCalibration/MultipleCameraPoseOptimization.h>
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
	const Ubitrack::Math::Vector< 2, T > m_maxDimension;
	
public:
	isPixelNotWithinScreen( const Ubitrack::Math::Vector< 2, T >& screenResolution )
		: m_maxDimension( screenResolution )
		{}
	
	bool operator()( Ubitrack::Math::Vector< 2, T > &pixel ) const
	{
		// just check if pixel is within given screen resolution.
		if ( pixel[ 0 ] < 0 || pixel [ 1 ] < 0 || pixel[ 0 ] >= m_maxDimension[ 0 ] || pixel [ 1 ] >= m_maxDimension[ 1 ] )
			return true;
		return false;
	}
};

}

template< typename T >
void TestMarkerBundleAdjustment( const std::size_t n_runs, const T epsilon )
{
	Vector< 2, T > screenResolution( 640, 480 );
	
	
	// run bundleAdjustment test for several times
	for( std::size_t i( 0 ); i< n_runs; ++i )
	{
		//change here for a really big bundle adjustment:
		const std::size_t n_p3d = Random::distribute_uniform( 10, 15 ); 
		const std::size_t n_cams = Random::distribute_uniform( 3, 20 );
		
		// random intrinsics matrix, never changes, assume always the same camera
		Matrix< 3, 3, T > cam( boost::numeric::ublas::identity_matrix< T >( 3 ) );
		cam( 0, 0 ) = Random::distribute_uniform< T >( 500, 800 );
		cam( 1, 1 ) = Random::distribute_uniform< T >( 500, 800 );
		cam( 0, 2 ) = screenResolution[ 0 ] * 0.5;
		cam( 1, 2 ) = screenResolution[ 1 ] * 0.5;
			
		typename Random::Vector< 2, T >::Normal randPixelNoise( 0, 0.5 ); // gaussian noise for 2d Pixels
		typename Random::Vector< 3, T >::Uniform randVector( -10.0, 10.0 ); // 3d points, assume meters as unit
		typename Random::Vector< 3, T >::Normal randPositionNoise( 0, 0.05 ); // gaussian noise for 3d points
		// typename Random::Vector< 3, T >::Uniform randPositionNoise( -0.25, 0.25 ); // uniform noise	for 3d points
		
		// random poses for the estimation
		typename Random::Quaternion< T >::Uniform randQuat;
		typename Random::Vector< 3, T >::Uniform randTranslation( -100, 100 ); //translation
		//alternatively one could use a random pose:
		// Random::Pose< T >::Uniform( -100, 100 );
		
	
		// generate 3d points (ground truth) + add some noise to points
		std::vector< Vector< 3, T > > points_3D;
		points_3D.reserve( n_p3d );
		std::generate_n ( std::back_inserter( points_3D ), n_p3d,  randVector );
		
		//generate some noise for the 3d points and add it
		std::vector< Vector< 3, T > > noisy_points_3D;
		noisy_points_3D.reserve( n_p3d );
		std::generate_n ( std::back_inserter( noisy_points_3D ), n_p3d,  randPositionNoise );
		std::transform( points_3D.begin(), points_3D.end(), noisy_points_3D.begin(), noisy_points_3D.begin(), std::plus< Vector< 3, T > >() );
		
		// std::cout << std::endl << std::endl;
		// std::copy( noisy_points_3D.begin(), noisy_points_3D.end(), std::ostream_iterator< Ubitrack::Math::Vector< 3, T > > ( std::cout, ", ") );
		// std::cout << std::endl;
		// std::copy( points_3D.begin(), points_3D.end(), std::ostream_iterator< Ubitrack::Math::Vector< 3, T > > ( std::cout, ", ") );
		
		// generate random poses and project on image plane
		std::vector< Pose > extrinsics;
		std::vector< Matrix< 3, 3, T > > intrinsics;
		std::vector< std::vector< Vector< 2, T > > > points_2D;
		std::vector < std::vector < Scalar < T > > > weights_2D;
		std::vector < Scalar < int > > localSizes;
		localSizes.reserve( n_cams );
		
		do
		{
			Quaternion rot( randQuat( ) );
			Vector< 3, T > trans ( randTranslation() );
			
			Pose pose6D( rot, trans );
			Matrix< 3, 4, T > proj( rot, trans );
			proj = boost::numeric::ublas::prod( cam, proj );
			
			std::vector< Vector< 2, T > > p2D;
			p2D.reserve( n_p3d );
			std::transform( points_3D.begin(), points_3D.end(), std::back_inserter( p2D ), Functors::ProjectVector< T >( proj ) );
			
			//generate some noise for the 2d Points
			std::vector< Vector< 2, T > > noisy_points_2D;
			noisy_points_2D.reserve( n_p3d );
			std::generate_n ( std::back_inserter( noisy_points_2D ), n_p3d,  randPixelNoise );
			std::transform( p2D.begin(), p2D.end(), noisy_points_2D.begin(), noisy_points_2D.begin(), std::plus< Vector< 2, T > >() );
			
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
				points_2D.push_back( noisy_points_2D );
				
				std::vector < Scalar < T > > weights;
				weights.reserve( n_2d );
				std::fill_n( std::back_inserter( weights ), n_2d, Scalar< T >( 1.0 ) );
				weights_2D.push_back( weights );
				
				localSizes.push_back( n_2d );
				//add some stuff here
			}
		}
		while( extrinsics.size() < n_cams );
	
		//prepare results
		const std::size_t minCorrespondences( 4 );
		std::vector < ErrorPose > errPoses;
		std::vector < Scalar< T > > poseWeights;
		
		// estimate result
		Ubitrack::Calibration::multipleCameraPoseEstimationWithLocalBundles ( points_3D
			, points_2D
			, weights_2D
			, extrinsics
			, intrinsics
			, minCorrespondences
			, errPoses
			, poseWeights
			, localSizes );
			
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
	//attention: will not work with floats so far.
	TestMarkerBundleAdjustment< double >( 10, 1e-3 );
}


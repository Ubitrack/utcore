

#include <utMath/Pose.h>
#include <utMath/Scalar.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Geometry/PointProjection.h>

#include <utMath/Random/Pose.h>
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../tools.h"

#include <utCalibration/BundleAdjustment.h>
//#include <utCalibration/MultipleCameraPoseOptimization.h> // <- no test yet

#include <iostream>
#include <algorithm>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>



using namespace Ubitrack::Math;

namespace { //anonymous namespace

template< typename T >
struct isPointWithinScreen
	: public std::binary_function< Ubitrack::Math::Matrix< T, 3, 4 >, Ubitrack::Math::Vector< T, 3 >, bool >
{
protected:
	const Ubitrack::Math::Vector< T, 2 > m_maxDimension;
	
public:
	isPointWithinScreen( const Ubitrack::Math::Vector< T, 2 >& screenResolution )
		: std::binary_function< Ubitrack::Math::Matrix< T, 3, 4 >, Ubitrack::Math::Vector< T, 3 >, bool >()
		, m_maxDimension( screenResolution )
		{}
	
	bool operator()( const Ubitrack::Math::Matrix< T, 3, 4 >& cam, Ubitrack::Math::Vector< T, 3 > &img2d ) const
	{
		const Ubitrack::Math::Vector< T, 2 > pixel = Ubitrack::Math::Geometry::ProjectPoint()( cam, img2d );
		// just check if pixel is within given screen resolution.
		if ( pixel[ 0 ] < 0 || pixel [ 1 ] < 0 || pixel[ 0 ] >= m_maxDimension[ 0 ] || pixel [ 1 ] >= m_maxDimension[ 1 ] )
			return false;
		return true;
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
		const std::size_t n_p3d = 30;//Random::distribute_uniform( 10, 15 ); 
		const std::size_t n_cams = Random::distribute_uniform( 3, 5 );
		
		// random intrinsics matrix, never changes, assume always the same camera
		Matrix< T, 3, 3 > cam( Matrix< T, 3, 3 >::identity() );
		cam( 0, 0 ) = Random::distribute_uniform< T >( 500, 800 );
		cam( 1, 1 ) = Random::distribute_uniform< T >( 500, 800 );
		cam( 0, 2 ) = screenResolution[ 0 ] * 0.5;
		cam( 1, 2 ) = screenResolution[ 1 ] * 0.5;
		cam( 2, 2 ) = 1;
			
		// typename Random::Vector< T, 2 >::Normal randPixelNoise( 0, 0.25 ); // gaussian noise for 2d Pixels -> usually very below 1 pixel
		typename Random::Vector< T, 3 >::Uniform randVector( -5.0, 5.0 ); // 3d points, assume meters as unit
		typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.05 ); // gaussian noise for 3d points -> can be quite high
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
		std::vector< Vector< T, 3 > > points_3D_noisy;
		points_3D_noisy.reserve( n_p3d );
		std::generate_n ( std::back_inserter( points_3D_noisy ), n_p3d,  randPositionNoise );
		std::transform( points_3D.begin(), points_3D.end(), points_3D_noisy.begin(), points_3D_noisy.begin(), std::plus< Vector< T, 3 > >() );
		
		// generate random poses and project on image plane
		std::vector< Pose > extrinsics_orig;
		std::vector< Pose > extrinsics_noisy;
		std::vector< std::vector< Vector< T, 2 > > > observed_points_2D;
		
		do
		{
			Quaternion rot( randQuat( ) );
			rot.normalize() ;
			
			Vector< T, 3 > trans ( randTranslation() );
			
			Pose pose6D( rot, trans );
			Matrix< T, 3, 4 > proj( rot, trans );
			
			std::vector< Vector< T, 2 > > points_2D;
			points_2D.reserve( n_p3d );
			Geometry::project_points( proj, points_3D.begin(), points_3D.end(), std::back_inserter( points_2D ) );
			
			//generate the 2D points without noise -> clear observation of ground truth data
			std::vector< Vector< T, 2 > > noisy_points_2D;
			noisy_points_2D.reserve( n_p3d );
			std::copy( points_2D.begin(), points_2D.end(), std::back_inserter( noisy_points_2D ) );
			
			
			// remove pixels outside the screen
			//noisy_points_2D.erase( std::remove_if( noisy_points_2D.begin(), noisy_points_2D.end(), isPixelNotWithinScreen< T >( screenResolution ) ), noisy_points_2D.end() ); 
			// const std::size_t n_2d( noisy_points_2D.size() );
			
			
			// at the moment count visible image points
			proj = boost::numeric::ublas::prod( cam, proj );
			const std::size_t n_2d = std::count_if( points_3D.begin(), points_3D.end(), std::bind1st( isPointWithinScreen< T >( screenResolution ), proj ) );
			
			if( n_2d > ( n_p3d /2 ) ) //enough points are visible?
			{
				// Pose noisyPose( rot, trans  );
				// Pose noisyPose( rot, trans + randPositionNoise() );
				const T rotEps( 0.01 );
				Pose noisyPose( Quaternion( rot.x() + Random::distribute_uniform< T >( -rotEps, rotEps )
						, rot.y() + Random::distribute_uniform< T >( -rotEps, rotEps )
						, rot.z() + Random::distribute_uniform< T >( -rotEps, rotEps )
						, rot.w() + Random::distribute_uniform< T >( -rotEps, rotEps ) )
						, trans + randPositionNoise() );
				
						
				extrinsics_orig.push_back( pose6D );
				extrinsics_noisy.push_back( noisyPose );
				observed_points_2D.push_back( noisy_points_2D );
			}
		}
		while( extrinsics_noisy.size() < n_cams );
		
		const T error3D = meanSummedDiff( points_3D_noisy, points_3D );
		const T errorPoseT = meanSummedTranslationDiff< T >( extrinsics_noisy, extrinsics_orig );
		const T errorPoseA = meanSummedAngularDiff< T >( extrinsics_noisy, extrinsics_orig );
		
		
		// std::cout << std::endl << std::endl << std::endl << std::endl;
		// std::copy( extrinsics_noisy.begin(), extrinsics_noisy.end(), std::ostream_iterator< Ubitrack::Math::Pose > ( std::cout, ", ") );	
		// std::cout << std::endl << std::endl;
		// std::copy( points_3D_noisy.begin(), points_3D_noisy.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 3 > > ( std::cout, ", ") );	
		// std::cout << std::endl << std::endl;
		
		Ubitrack::Calibration::simpleBundleAdjustment( observed_points_2D, extrinsics_noisy, points_3D_noisy );
		// call to templated function (does not link on windows):
		// Ubitrack::Calibration::simpleBundleAdjustment( observed_points_2D.begin(), observed_points_2D.end(), intrinsics.begin(), extrinsics.begin(), points_3D_noisy.begin(), points_3D_noisy.end() );
		
		
		const T optError3D = meanSummedDiff( points_3D_noisy, points_3D );
		const T optErrorPoseT = meanSummedTranslationDiff< T >( extrinsics_noisy, extrinsics_orig );
		const T optErrorPoseA = meanSummedAngularDiff< T >( extrinsics_noisy, extrinsics_orig );
		
		
		
		BOOST_WARN_MESSAGE( error3D >= optError3D, "Error 3D point        " << error3D << " vs. " << optError3D << " (before vs. later)" );
		BOOST_WARN_MESSAGE( errorPoseT >= optErrorPoseT, "Error (pose) position " << errorPoseT << " vs. " << optErrorPoseT << " (before vs. later)" );
		BOOST_WARN_MESSAGE( errorPoseA >= optErrorPoseA, "Error (pose) angle    " << errorPoseA << " vs. " << optErrorPoseA << " (before vs. later)" );
		
		
		// print out 3d values
		// std::copy( points_3D_noisy.begin(), points_3D_noisy.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 3 > > ( std::cout, ", ") );	
		// std::cout << std::endl << std::endl;
		// std::copy( points_3D.begin(), points_3D.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 3 > > ( std::cout, ", ") );	
		// std::cout << std::endl << std::endl;
		
		// print out poses
		// std::copy( extrinsics_noisy.begin(), extrinsics_noisy.end(), std::ostream_iterator< Ubitrack::Math::Pose > ( std::cout, ", ") );	
		// std::cout << std::endl << std::endl;
		// std::copy( extrinsics_orig.begin(), extrinsics_orig.end(), std::ostream_iterator< Ubitrack::Math::Pose > ( std::cout, ", ") );	
		// std::cout << std::endl << std::endl;
	}	
};

void TestBundleAdjustment()
{
	// attention: works also with float now :)
	TestMarkerBundleAdjustment< double >( 10, 1e-3 );
	TestMarkerBundleAdjustment< float >( 10, 1e-3 );
}


#include <utAlgorithm/3DPointReconstruction.h>
#include <utMath/Geometry/PointProjection.h>
#include <utMath/Stochastic/identity_iterator.h>

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include "../tools.h"

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack;

template< typename T >
void Test2Cameras( const std::size_t n_runs, const T epsilon )
{
	typename Math::Random::Quaternion< T >::Uniform randQuat;
	typename Math::Random::Vector< T, 3 >::Uniform randTranslation( -10., 10. ); //translation
	
	typename Math::Random::Vector< T, 3 >::Uniform randVector( -1., 1. ); // 3d Points
	
	for( std::size_t j=0; j<n_runs; ++j )
	{

		Math::Pose CamPose1( randQuat() , randTranslation() );
		Math::Pose CamPose2( randQuat() , randTranslation() );
		Math::Matrix< T, 3, 4 > proj1( CamPose1 );
		Math::Matrix< T, 3, 4 > proj2( CamPose2 );

		Math::Matrix< T, 3, 3 > K = Math::Matrix< T, 3, 3 >::identity();
		K( 0, 0 ) = K( 1, 1 ) = 500;
		K( 0, 2 ) = 320;
		K( 1, 2 ) = 240;

		// generate the projections
		proj1 = boost::numeric::ublas::prod( K, proj1 );
		proj2 = boost::numeric::ublas::prod( K, proj2 );
		
		//compute random points
		const std::size_t n( Math::Random::distribute_uniform< std::size_t >( 10, 30 ) );
		std::vector< Math::Vector< T, 3 > > objPoints;
		objPoints.reserve( n );
		std::generate_n ( std::back_inserter( objPoints ), n,  randVector );
		
		//project the points onto the first image screen
		std::vector< Math::Vector< T, 2 > > points1;
		points1.reserve( n );
		Math::Geometry::project_points( proj1, objPoints.begin(), objPoints.end(), std::back_inserter( points1 ) );
		
		//project the points onto the second image screen
		std::vector< Math::Vector< T, 2 > > points2;
		points2.reserve( n );
		Math::Geometry::project_points( proj2, objPoints.begin(), objPoints.end(), std::back_inserter( points2 ) );
		
		for( std::size_t i=0; i<n; ++i )
		{
			//check now the reconstructed points		
			Math::Vector< T, 3 > p3D = Algorithm::get3DPosition( proj1, proj2, points1[ i ], points2[ i ] );
			const T diffError = vectorDiff( p3D, objPoints[ i ] );
			BOOST_CHECK( diffError < epsilon );	
		}
	}
}


template< typename T >
void TestMulitpleCameras( const std::size_t n_runs , const T epsilon )
{
	typename Math::Random::Quaternion< T >::Uniform randQuat;
	typename Math::Random::Vector< T, 3 >::Uniform randTranslation( -10., 10. ); //translation
	
	typename Math::Random::Vector< T, 3 >::Uniform randVector( -1., 1. ); // 3d Points
	// typename Random::Vector< T, 3 >::Normal randPositionNoise( 0, 0.2 ); // translation gaussian noise
	
	// check for multi-camera setup
	for( std::size_t j=0; j<n_runs; ++j )
	{
		// determine the amount of cameras to be used
		const std::size_t n_cams( Math::Random::distribute_uniform< std::size_t >( 2, 10 ) );
		
		// intialize some cameras:
		std::vector< Math::Matrix< T, 3, 4 > > matrices; // stores projection matrices
		matrices.reserve( n_cams );
		
		for( std::size_t i( 0 ); i < n_cams; ++i )
		{
			Math::Pose CamPose( randQuat() , randTranslation() );
			Math::Matrix< T, 3, 4 > p( CamPose );
			matrices.push_back( p );
		}			
			
		// generate a random vector
		const Math::Vector< T, 3 >  randVec = randVector();
		
		// generate some 2D observation of 3D random vector
		std::vector< Math::Vector< T, 2 > > pts; //stores image points
		pts.reserve( n_cams );
		std::transform( matrices.begin(), matrices.end(), Util::identity< Math::Vector< T, 3 > >( randVec ).begin()
			, std::back_inserter( pts ), Math::Geometry::ProjectPoint() );

		//estimate final result
		Math::Vector< T, 3 > p3D = Algorithm::get3DPosition( matrices, pts, 0 );
		const T diffError = vectorDiff( p3D, randVec );
		BOOST_CHECK_SMALL( diffError, epsilon );
	}
}

void Test3DPointReconstruction()
{
	Test2Cameras< float >( 1000, 1e-2f );
	Test2Cameras< double >( 1000, 1e-3 );
	TestMulitpleCameras< float >( 1000, 1e-2f );
	TestMulitpleCameras< double >( 1000, 1e-3 );
}
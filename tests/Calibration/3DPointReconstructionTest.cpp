#include <utCalibration/3DPointReconstruction.h>
#include <utMath/Functors/Vector3Functors.h>
#include "../tools.h"

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack;
namespace ublas = boost::numeric::ublas;

void Test3DPointReconstruction()
{
	for( int j=0; j<100; ++j )
	{

		Math::Pose CamPose1( randomQuaternion() , randomVector< 3, double >() );
		Math::Pose CamPose2( randomQuaternion() , randomVector< 3, double >() );
		Math::Matrix< 3, 4, double > p1( CamPose1 );
		Math::Matrix< 3, 4, double > p2( CamPose2 );

		Math::Matrix< 3, 3, double > K = Math::Matrix< 3, 3, double >::identity();
		K( 0, 0 ) = K( 1, 1 ) = 500;
		K( 0, 2 ) = 320;
		K( 1, 2 ) = 240;

		p1 = ublas::prod( K, p1 );
		p2 = ublas::prod( K, p2 );
		
		//compute random points
		unsigned int num_correspondences = static_cast< unsigned > ( random( 10.0, 30.0 ) );
		std::vector< Math::Vector< 3 > > objPoints;
		std::vector< Math::Vector< 2 > > points1;
		std::vector< Math::Vector< 2 > > points2;
		objPoints.reserve( num_correspondences );
		points1.reserve( num_correspondences );
		points2.reserve( num_correspondences );
		
		//project the points onto the image screen
		for( unsigned int i=0; i<num_correspondences; ++i )
		{
			Math::Vector< 3 > v = randomVector< 3, double >();
			Math::Vector< 2 > v1 = Math::Functors::project3x4_vector3< double >()( p1, v );
			Math::Vector< 2 > v2 = Math::Functors::project3x4_vector3< double >()( p2, v );

			objPoints.push_back( v );
			points1.push_back( Math::Vector< 2 >( v1 ) );
			points2.push_back( Math::Vector< 2 >( v2 ) );
		
			//check now the reconstructed points		
			Math::Vector< 3 > tmp = Calibration::get3DPosition( p1, p2, points1[ i ], points2[ i ] );
			BOOST_CHECK_SMALL( vectorDiff( tmp, objPoints[ i ] ), 1e-05 );	
		}
	}
	
	// check for multi-camera setup
	for( int j=0; j<100; ++j )
	{

		unsigned int num_cameras = static_cast< unsigned > ( random( 2.0, 10.0 ) );
		std::vector< Math::Matrix< 3, 4, double > > matrices; // stores projection matrices
		std::vector< Math::Vector< 2 , double > > pts; //stores image points
		matrices.reserve( num_cameras );
		pts.reserve( num_cameras );

		// random vector
		Math::Vector< 3 > v = randomVector< 3, double >() * random( 2.0, 100.0 );

		for( unsigned i( 0 ); i < num_cameras; ++i )
		{
			Math::Pose CamPose( randomQuaternion() , randomVector< 3, double >() );
			Math::Matrix< 3, 4, double > p( CamPose );
			Math::Vector< 2 > p2d = Math::Functors::project3x4_vector3< double >()( p, v );
			
			matrices.push_back( p );
			pts.push_back( p2d );
		}
		
		Math::Vector< 3, double > p3D = Calibration::get3DPosition( matrices, pts, 0 );
		BOOST_CHECK_SMALL( vectorDiff( p3D, v ), 1e-05 );
	}
	
}
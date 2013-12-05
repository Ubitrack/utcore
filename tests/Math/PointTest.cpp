

#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>

#include <utMath/Geometry/PointProjection.h>
#include <utMath/Geometry/PointTransformation.h>

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack::Math;

//right now this function does only compilation test, nothing more 
template< typename T >
void testBasicPointTransformations( const std::size_t n )
{
	
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< 3, T >::Uniform randTranslation( -10, 10 ); //translation

	// random pose
	Quaternion rot( randQuat( ) );
	Vector< 3, T > trans ( randTranslation() );
	
	// transformation matrices
	Matrix< 4, 4, T > mat4x4( rot, trans );
	Matrix< 3, 4, T > mat3x4( rot, trans );
	Matrix< 3, 3, T > mat3x3( rot );
	// following matrix set to identity -> boring
	Matrix< 2, 3, T > mat2x3; 
	mat2x3( 0, 0 ) = mat2x3( 1, 1 ) =  1;
	mat2x3( 0, 1 ) = mat2x3( 1, 0 ) =  0;
	mat2x3( 0, 2 ) = trans( 0 );
	mat2x3( 1, 2 ) = trans( 1 );
	mat2x3( 2, 2 ) = trans( 2 );
	{
		// random 3D parameters (as 2D vector, last dimension as one)
		typename Random::Vector< 2, T >::Uniform randPoints2D( -5, 5 ); 

		// generate some random points
		std::vector< Ubitrack::Math::Vector< 2, T > > pts3D2;
		pts3D2.reserve( n );
		std::generate_n ( std::back_inserter( pts3D2 ), n,  randPoints2D );

		std::vector< Ubitrack::Math::Vector< 2, T > > pts3DOut11;
		pts3DOut11.reserve( n );
		Geometry::transform_points( mat2x3, pts3D2.begin(), pts3D2.end(), std::back_inserter( pts3DOut11 ) );
		
		std::vector< Ubitrack::Math::Vector< 3, T > > pts3DOut12;
		pts3DOut12.reserve( n );
		Geometry::transform_points( mat3x3, pts3D2.begin(), pts3D2.end(), std::back_inserter( pts3DOut12 ) );
		
		std::vector< Ubitrack::Math::Vector< 3, T > > pts3DOut13;
		pts3DOut13.reserve( n );
		Geometry::transform_points( mat3x4, pts3D2.begin(), pts3D2.end(), std::back_inserter( pts3DOut13 ) );
		
		std::vector< Ubitrack::Math::Vector< 4, T > > pts3DOut14;
		pts3DOut14.reserve( n );
		Geometry::transform_points( mat4x4, pts3D2.begin(), pts3D2.end(), std::back_inserter( pts3DOut14 ) );
		
	}
	
	
	
	{
		// random 3D parameters (as 3D vector, last dimension as one)
		typename Random::Vector< 3, T >::Uniform randPoints3D( -5, 5 ); 

		// generate some random points
		std::vector< Ubitrack::Math::Vector< 3, T > > pts3D1;
		pts3D1.reserve( n );
		std::generate_n ( std::back_inserter( pts3D1 ), n,  randPoints3D );

		std::vector< Ubitrack::Math::Vector< 2, T > > pts3DOut11;
		pts3DOut11.reserve( n );
		Geometry::transform_points( mat2x3, pts3D1.begin(), pts3D1.end(), std::back_inserter( pts3DOut11 ) );
		
		std::vector< Ubitrack::Math::Vector< 3, T > > pts3DOut12;
		pts3DOut12.reserve( n );
		Geometry::transform_points( mat3x3, pts3D1.begin(), pts3D1.end(), std::back_inserter( pts3DOut12 ) );
		
		std::vector< Ubitrack::Math::Vector< 3, T > > pts3DOut13;
		pts3DOut13.reserve( n );
		Geometry::transform_points( mat3x4, pts3D1.begin(), pts3D1.end(), std::back_inserter( pts3DOut13 ) );
		
		std::vector< Ubitrack::Math::Vector< 4, T > > pts3DOut14;
		pts3DOut14.reserve( n );
		Geometry::transform_points( mat4x4, pts3D1.begin(), pts3D1.end(), std::back_inserter( pts3DOut14 ) );
		
	}
	
	
	{
		// random 3D parameters (as 4D vector)
		typename Random::Vector< 4, T >::Uniform randPoints4D( -5, 5 ); 

		// generate some random points
		std::vector< Ubitrack::Math::Vector< 4, T > > pts3D2;
		pts3D2.reserve( n );
		std::generate_n ( std::back_inserter( pts3D2 ), n,  randPoints4D );

		std::vector< Ubitrack::Math::Vector< 3, T > > pts3DOut11;
		pts3DOut11.reserve( n );
		
		std::vector< Ubitrack::Math::Vector< 3, T > > pts3DOut13;
		pts3DOut13.reserve( n );
		Geometry::transform_points( mat3x4, pts3D2.begin(), pts3D2.end(), std::back_inserter( pts3DOut13 ) );
		
		std::vector< Ubitrack::Math::Vector< 4, T > > pts3DOut14;
		pts3DOut14.reserve( n );
		Geometry::transform_points( mat4x4, pts3D2.begin(), pts3D2.end(), std::back_inserter( pts3DOut14 ) );
		
	}
}

//right now this function does only compilation test, nothing more as well
template< typename T >
void testBasicPointProjection( const std::size_t n )
{
	const Vector< 2, T > screenResolution( 640, 480 );
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< 3, T >::Uniform randTranslation( -10, 10 ); //translation

	// random pose
	Quaternion rot( randQuat( ) );
	Vector< 3, T > trans ( randTranslation() );
	
	// random intrinsics matrix, never changes, assume always the same camera
	Matrix< 3, 3, T > cam( Matrix< 3, 3, T >::identity() );
	cam( 0, 0 ) = Random::distribute_uniform< T >( 500, 800 );
	cam( 1, 1 ) = Random::distribute_uniform< T >( 500, 800 );
	// take care of ubitracks camera interpretation -> last column negative entries
	cam( 0, 2 ) = -screenResolution[ 0 ] * 0.5;
	cam( 1, 2 ) = -screenResolution[ 1 ] * 0.5;
	cam( 2, 2 ) = -1;
	
	// projection to 2D image plane
	Matrix< 3, 4, T > projection( rot, trans );
	projection = boost::numeric::ublas::prod( cam, projection );
	
	
	
	{
		// random 3D parameters (as 2D vector)
		typename Random::Vector< 2, T >::Uniform randPoints2D( -5, 5 ); 

		// generate some random points
		std::vector< Ubitrack::Math::Vector< 2, T > > pts3D2;
		pts3D2.reserve( n );
		std::generate_n ( std::back_inserter( pts3D2 ), n,  randPoints2D );

		std::vector< Ubitrack::Math::Vector< 2, T > > pts2DOut1;
		pts2DOut1.reserve( n );
		Geometry::project_points( projection, pts3D2.begin(), pts3D2.end(), std::back_inserter( pts2DOut1 ) );
	}
	
	
	{
		// random 3D parameters (as 3D vector)
		typename Random::Vector< 3, T >::Uniform randPoints3D( -5, 5 ); 

		// generate some random points
		std::vector< Ubitrack::Math::Vector< 3, T > > pts3D3;
		pts3D3.reserve( n );
		std::generate_n ( std::back_inserter( pts3D3 ), n,  randPoints3D );
		
		std::vector< Ubitrack::Math::Vector< 2, T > > pts2DOut2;
		pts2DOut2.reserve( n );
		Geometry::project_points( projection, pts3D3.begin(), pts3D3.end(), std::back_inserter( pts2DOut2 ) );
	}
	
	{
		// random 3D parameters (as 3D vector)
		typename Random::Vector< 4, T >::Uniform randPoints4D( -5, 5 ); 

		// generate some random points
		std::vector< Ubitrack::Math::Vector< 4, T > > pts3D4;
		pts3D4.reserve( n );
		std::generate_n ( std::back_inserter( pts3D4 ), n,  randPoints4D );
		
		std::vector< Ubitrack::Math::Vector< 2, T > > pts2DOut3;
		pts2DOut3.reserve( n );
		Geometry::project_points( projection, pts3D4.begin(), pts3D4.end(), std::back_inserter( pts2DOut3 ) );
	}
	
	//just to get rid of compilation warning
	BOOST_CHECK_SMALL( 0.01, 0.02 );
}

void TestPoints()
{
	testBasicPointTransformations< float >( 10000 );
	testBasicPointTransformations< double >( 10000 );
	testBasicPointProjection< float >( 10000 );
	testBasicPointProjection< double >( 10000 );
}



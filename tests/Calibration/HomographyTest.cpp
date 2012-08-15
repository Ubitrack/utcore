#include <boost/test/unit_test.hpp>
#include <math.h>
#include <boost/test/floating_point_comparison.hpp>
#include <utCalibration/Homography.h>
#include "../tools.h"
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <iostream>

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

void TestHomographyDLT()
{
	// check if standard corners return identity
	std::vector< Vector< 2, float > > stdCorners( 4 );
	for ( int i = 0; i < 4; i++ )
	{
		stdCorners[ i ][ 0 ] = ( i & 2 )         ? 0.5f : -0.5f;
		stdCorners[ i ][ 1 ] = ( ( i + 1 ) & 2 ) ? -0.5f : 0.5f;
	}
	
	Matrix< 3, 3, float > H( Ubitrack::Calibration::homographyDLT( stdCorners, stdCorners ) );
	BOOST_CHECK_SMALL( homMatrixDiff( H, ublas::identity_matrix< float >( 3 ) ), 0.001 );

	// do some tests with random homographies
	for ( int iTest = 0; iTest < 100; iTest++ )
	{
		// create a random homography
		Matrix< 3, 3, float > Htest;
		randomMatrix( Htest );
				
		// create & transform random points
		// use at least 10 correspondences, as randomness may lead to poorly conditioned problems...
		unsigned pointNum = 10 + ( iTest % 10 ); 
		std::vector< Vector< 2, float > > fromPoints( pointNum );
		std::vector< Vector< 2, float > > toPoints( pointNum );
		for ( unsigned i = 0; i < pointNum; i++ )
		{
			Vector< 3, float > x;
			x( 0 ) = fromPoints[ i ]( 0 ) = random( -100.0f, 100.0f );
			x( 1 ) = fromPoints[ i ]( 1 ) = random( -100.0f, 100.0f );
			x( 2 ) = 1.0f;
			
			Vector< 3, float > xp = ublas::prod( Htest, x );
			toPoints[ i ] = ublas::subrange( xp, 0, 2 ) / xp( 2 );
		}
		
		H = Ubitrack::Calibration::homographyDLT( fromPoints, toPoints );

		BOOST_CHECK_SMALL( homMatrixDiff( H, Htest ), 1e-2 );
	}
}


void TestSquareHomography()
{
	// check if standard corners return identity
	std::vector< Vector< 2, float > > stdCorners( 4 );
	for ( int i = 0; i < 4; i++ )
	{
		stdCorners[ i ][ 0 ] = ( i & 2 )         ? 0.5f : -0.5f;
		stdCorners[ i ][ 1 ] = ( ( i + 1 ) & 2 ) ? -0.5f : 0.5f;
	}
	
	Matrix< 3, 3, float > H( Ubitrack::Calibration::squareHomography( stdCorners ) );
	BOOST_CHECK_SMALL( homMatrixDiff( H, ublas::identity_matrix< float >( 3 ) ), 1e-3 );
	
	// do some tests with random homographies
	for ( int iTest = 0; iTest < 100; iTest++ )
	{
		// create a random homography
		Matrix< 3, 3, float > Htest;
		randomMatrix( Htest );
				
		// transform standard corners by Htest
		std::vector< Vector< 2, float > > tCorners( 4 );
		for ( int i = 0; i < 4; i++ )
		{
			Vector< 3, float > x;
			ublas::subrange( x, 0, 2 ) = stdCorners[ i ];
			x( 2 ) = 1.0f;
			
			Vector< 3, float > xp = ublas::prod( Htest, x );
			tCorners[ i ] = ublas::subrange( xp, 0, 2 ) / xp( 2 );
		}
		
		H = Ubitrack::Calibration::squareHomography( tCorners );

		BOOST_CHECK_SMALL( homMatrixDiff( H, Htest ), 1e-2 );
	}
}


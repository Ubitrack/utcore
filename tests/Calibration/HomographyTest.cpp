
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Random/Vector.h>
#include <utCalibration/Homography.h>

#include "../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;

template< typename T >
void TestHomographyDLTIdentity( const T epsilon )
{
	// check if standard corners return identity
	std::vector< Vector< 2, T > > stdCorners( 4 );
	for ( std::size_t i = 0; i < 4; i++ )
	{
		stdCorners[ i ][ 0 ] = ( i & 2 )         ? 0.5 : -0.5;
		stdCorners[ i ][ 1 ] = ( ( i + 1 ) & 2 ) ? -0.5 : 0.5;
	}
	
	Matrix< 3, 3, T > H( Ubitrack::Calibration::homographyDLT( stdCorners, stdCorners ) );
	BOOST_CHECK_SMALL( homMatrixDiff( H, Matrix< 3, 3, T >::identity() ), epsilon );
}

template< typename T >
void TestSquareHomography( const std::size_t n_runs, const T epsilon )
{
	// check if standard corners return identity
	std::vector< Vector< 2, T > > stdCorners( 4 );
	for ( std::size_t i = 0; i < 4; i++ )
	{
		stdCorners[ i ][ 0 ] = ( i & 2 )         ? 0.5f : -0.5f;
		stdCorners[ i ][ 1 ] = ( ( i + 1 ) & 2 ) ? -0.5f : 0.5f;
	}
	
	Matrix< 3, 3, T > H( Ubitrack::Calibration::squareHomography( stdCorners ) );
	BOOST_CHECK_SMALL( homMatrixDiff( H, Matrix< 3, 3, T >::identity() ), epsilon );
	
	
	for ( std::size_t iTest = 0; iTest < n_runs; iTest++ )
	{
		Matrix< 3, 3, T > Htest;
		randomMatrix( Htest );
					
		// transform standard corners by Htest
		std::vector< Vector< 2, T > > tCorners( 4 );
		for ( std::size_t i = 0; i < 4; i++ )
		{
			Vector< 3, T > x( stdCorners[ i ] ( 0 ), stdCorners[ i ] ( 1 ), 1. );
			Vector< 3, T > xp = ublas::prod( Htest, x );
			tCorners[ i ] = ublas::subrange( xp, 0, 2 ) / xp( 2 );
		}
			
		H = Ubitrack::Calibration::squareHomography( tCorners );
		BOOST_CHECK_SMALL( homMatrixDiff( H, Htest ), epsilon );
	}
}
template < typename T >
void TestHomographyDLT( const std::size_t n_runs, const T epsilon )	
{
	// initialize random number generation.
	typename Random::Vector< 2, T >::Uniform randVector( -100, 100 ); // some 3D points
	
	// do some tests with random homographies
	for ( std::size_t iTest = 0; iTest < n_runs; iTest++ )
	{
		// create a random homography
		Matrix< 3, 3, T > Htest;
		randomMatrix( Htest );
				
		// create & transform random points
		// use at least 10 correspondences, as randomness may lead to poorly conditioned problems...
		const std::size_t n( Random::distribute_uniform< std::size_t >( 10, 50 ) );
		
		std::vector< Vector< 2, T > > fromPoints;
		fromPoints.reserve( n );
		std::generate_n ( std::back_inserter( fromPoints ), n,  randVector );
		
		std::vector< Vector< 2, T > > toPoints( n );
	
		for ( std::size_t i = 0; i < n; ++i )
		{
			Vector< 3, T > x( fromPoints[ i ]( 0 ), fromPoints[ i ]( 1 ), 1. );
			Vector< 3, T > xp = ublas::prod( Htest, x );
			toPoints[ i ] = ublas::subrange( xp, 0, 2 ) / xp( 2 );
		}
		
		Matrix< 3, 3, T > H = Ubitrack::Calibration::homographyDLT( fromPoints, toPoints );

		BOOST_CHECK_SMALL( homMatrixDiff( H, Htest ), epsilon );
	}
}

void TestHomography()
{
	TestHomographyDLTIdentity< double >( 1e-6 );
	TestSquareHomography< double >( 1000, 1e-6 );
	TestHomographyDLT< double >( 1000, 1e-6 );
	
	TestHomographyDLTIdentity< float >( 1e-3 );
	TestSquareHomography< float >( 1000, 1e-2 );
	TestHomographyDLT< float >( 1000, 1e-2 );
}

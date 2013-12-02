#include <boost/test/unit_test.hpp>
#include <math.h>
#include <boost/test/floating_point_comparison.hpp>
#include <utCalibration/Projection.h>
#include "../tools.h"
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <iostream>

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

Matrix< 3, 4, float > randomProjection( Matrix< 3, 3, float >& k, Matrix< 3, 3, float >& r, Vector< 3, float >& t )
{
	// create a random rotation
	Quaternion q( random( -1.0, 1.0 ), random( -1.0, 1.0 ), random( -1.0, 1.0 ), random( -1.0, 1.0 ) );
	q.normalize();
	r = Matrix< 3, 3, float >( q );

	// create random translation
	t( 0 ) = random( -1000.0f, 1000.0f );
	t( 1 ) = random( -1000.0f, 1000.0f );
	t( 2 ) = random( -1000.0f, 0.01f );

	// create rt matrix
	Matrix< 3, 4, float > rt;
	ublas::subrange( rt, 0, 3, 0, 3 ) = r;
	ublas::column( rt, 3 ) = t;

	// create random intrinsics matrix
	k( 1, 0 ) = k( 2, 0 ) = k( 2, 1 ) = 0.0f; k( 2, 2 ) = -1.0f;
	k( 0, 0 ) = random( 1.0f, 1000.0f ); k( 1, 1 ) = random( 1.0f, 1000.0f );
	k( 0, 1 ) = random( -200.0f, 200.0f );
	k( 0, 2 ) = random( -200.0f, 200.0f ); k( 1, 2 ) = random( -200.0f, 200.0f );

	// multiply k rt and return
	return ublas::prod( k, rt );
}


Matrix< 3, 4, float > randomProjection()
{ Matrix< 3, 3, float > k; Matrix< 3, 3, float > r; Vector< 3, float > t; return randomProjection( k, r, t ); }


void TestProjectionDLT()
{
	// do some tests with random projections
	for ( int iTest = 0; iTest < 100; iTest++ )
	{
		// create a random projection matrix
		Matrix< 3, 4, float > Ptest( randomProjection() );

		// create & transform random points
		// use at least 20 correspondences, as randomness may lead to poorly conditioned problems...
		unsigned pointNum = 20 + ( iTest % 10 );
		std::vector< Vector< 3, float > > fromPoints( pointNum );
		std::vector< Vector< 2, float > > toPoints( pointNum );
		for ( unsigned i = 0; i < pointNum; i++ )
		{
			Vector< 4, float > x;
			x( 0 ) = fromPoints[ i ]( 0 ) = random( -100.0f, 100.0f );
			x( 1 ) = fromPoints[ i ]( 1 ) = random( -100.0f, 100.0f );
			x( 2 ) = fromPoints[ i ]( 2 ) = random( -100.0f, 100.0f );
			x( 3 ) = 1.0f;

			Vector< 3, float > xp = ublas::prod( Ptest, x );
			toPoints[ i ] = ublas::subrange( xp, 0, 2 ) / xp( 2 );
		}

		Matrix< 3, 4, float > P( Ubitrack::Calibration::projectionDLT( fromPoints, toPoints ) );

		BOOST_CHECK_SMALL( homMatrixDiff( P, Ptest ), 1e-3f );
	}
}


void TestDecomposeProjection()
{
	// do some tests with random projections
	for ( int iTest = 0; iTest < 100; iTest++ )
	{
		// create random projection
		Matrix< 3, 3, float > k;
		Matrix< 3, 3, float > r;
		Vector< 3, float > t;
		Matrix< 3, 4, float > p = randomProjection( k, r, t );

		// decompose
		Matrix< 3, 3, float > kEst;
		Matrix< 3, 3, float > rEst;
		Vector< 3, float > tEst;
		Ubitrack::Calibration::decomposeProjection( kEst, rEst, tEst, p );
		
		// check
		BOOST_CHECK_SMALL( matrixDiff( k, kEst ), 1e-3 );
		BOOST_CHECK_SMALL( matrixDiff( r, rEst ), 1e-3 );
		BOOST_CHECK_SMALL( vectorDiff( t, tEst ), 1e-3 );
	}
}


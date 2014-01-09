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

Matrix< float, 3, 4 > randomProjection( Matrix< float, 3, 3 >& k, Matrix< float, 3, 3 >& r, Vector< float, 3 >& t )
{
	// create a random rotation
	Quaternion q( random( -1.0, 1.0 ), random( -1.0, 1.0 ), random( -1.0, 1.0 ), random( -1.0, 1.0 ) );
	q.normalize();
	r = Matrix< float, 3, 3 >( q );

	// create random translation
	t( 0 ) = random( -1000.0f, 1000.0f );
	t( 1 ) = random( -1000.0f, 1000.0f );
	t( 2 ) = random( -1000.0f, 0.01f );

	// create rt matrix
	Matrix< float, 3, 4 > rt;
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


Matrix< float, 3, 4 > randomProjection()
{ Matrix< float, 3, 3 > k; Matrix< float, 3, 3 > r; Vector< float, 3 > t; return randomProjection( k, r, t ); }


void TestProjectionDLT()
{
	// do some tests with random projections
	for ( int iTest = 0; iTest < 100; iTest++ )
	{
		// create a random projection matrix
		Matrix< float, 3, 4 > Ptest( randomProjection() );

		// create & transform random points
		// use at least 20 correspondences, as randomness may lead to poorly conditioned problems...
		unsigned pointNum = 20 + ( iTest % 10 );
		std::vector< Vector< float, 3 > > fromPoints( pointNum );
		std::vector< Vector< float, 2 > > toPoints( pointNum );
		for ( unsigned i = 0; i < pointNum; i++ )
		{
			Vector< float, 4 > x;
			x( 0 ) = fromPoints[ i ]( 0 ) = random( -100.0f, 100.0f );
			x( 1 ) = fromPoints[ i ]( 1 ) = random( -100.0f, 100.0f );
			x( 2 ) = fromPoints[ i ]( 2 ) = random( -100.0f, 100.0f );
			x( 3 ) = 1.0f;

			Vector< float, 3 > xp = ublas::prod( Ptest, x );
			toPoints[ i ] = ublas::subrange( xp, 0, 2 ) / xp( 2 );
		}

		Matrix< float, 3, 4 > P( Ubitrack::Calibration::projectionDLT( fromPoints, toPoints ) );

		BOOST_CHECK_SMALL( homMatrixDiff( P, Ptest ), 1e-3f );
	}
}


void TestDecomposeProjection()
{
	// do some tests with random projections
	for ( int iTest = 0; iTest < 100; iTest++ )
	{
		// create random projection
		Matrix< float, 3, 3 > k;
		Matrix< float, 3, 3 > r;
		Vector< float, 3 > t;
		Matrix< float, 3, 4 > p = randomProjection( k, r, t );

		// decompose
		Matrix< float, 3, 3 > kEst;
		Matrix< float, 3, 3 > rEst;
		Vector< float, 3 > tEst;
		Ubitrack::Calibration::decomposeProjection( kEst, rEst, tEst, p );
		
		// check
		BOOST_CHECK_SMALL( matrixDiff( k, kEst ), 1e-3f );
		BOOST_CHECK_SMALL( matrixDiff( r, rEst ), 1e-3f );
		BOOST_CHECK_SMALL( vectorDiff( t, tEst ), 1e-3f );
	}
}


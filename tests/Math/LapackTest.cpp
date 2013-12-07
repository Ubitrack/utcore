#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <utMath/Matrix.h>

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/lapack/gels.hpp>
#endif

using namespace Ubitrack::Math;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

#ifdef HAVE_LAPACK
void TestBasicEigenvalues()
{
	using namespace Ubitrack;

	Math::Matrix< double, 3, 3 > A;
	A ( 0, 0 ) = 1; A ( 0, 1 ) = 2; A ( 0, 2 ) = 3;
	                A ( 1, 1 ) = 2; A ( 1, 2 ) = 4;
					                A ( 2, 2 ) = 5;

	Math::Vector< double, 3 > W;
	lapack::syev ( 'V', 'U', A, W, lapack::optimal_workspace() );

	Math::Matrix< double, 3, 3 > Atruth;
	Math::Vector< double, 3 > Wtruth;

	Atruth ( 0, 0 ) =  0.23319198; Atruth ( 0, 1 ) =  0.88765034; Atruth ( 0, 2 ) =  0.39711255;
	Atruth ( 1, 0 ) =  0.73923874; Atruth ( 1, 1 ) = -0.42713229; Atruth ( 1, 2 ) =  0.52065737;
	Atruth ( 2, 0 ) = -0.63178128; Atruth ( 2, 1 ) = -0.17214786; Atruth ( 2, 2 ) =  0.75578934;

	Wtruth[ 0 ] = -0.78765311;
	Wtruth[ 1 ] = -0.54419831;
	Wtruth[ 2 ] =  9.33185141;

	for ( int i = 0; i < 3; ++i )
	{
		for ( int j = 0; j < 3; ++j )
		{
			BOOST_CHECK_CLOSE ( A ( i, j ), Atruth ( i, j ), 10e-6 );
		}
		BOOST_CHECK_CLOSE ( W ( i ), Wtruth ( i ), 10e-6 );
	}
}

void TestLinearLeastSquares()
{
	using namespace Ubitrack;

	Math::Matrix< double, 10, 2 > A;
	Math::Matrix< double, 10, 1 > b;
	for (int i = 0; i < 10; ++i)
	{
		A( i, 0 ) = 1;
		A( i, 1 ) = i+1;
	}

	b( 0, 0 ) = 444;
	b( 1, 0 ) = 458;
	b( 2, 0 ) = 478;
	b( 3, 0 ) = 493;
	b( 4, 0 ) = 506;
	b( 5, 0 ) = 516;
	b( 6, 0 ) = 523;
	b( 7, 0 ) = 531;
	b( 8, 0 ) = 543;
	b( 9, 0 ) = 571;

	lapack::gels( 'N', A, b );

 	BOOST_CHECK_CLOSE( b( 0, 0 ), 436.2, 10e-6 );
 	BOOST_CHECK_CLOSE( b( 1, 0 ), 12.7454545454545, 10e-6 );
}


#endif // HAVE_LAPACK

void TestLapack()
{
#ifdef HAVE_LAPACK
	// check if singular value decomposition works
	Matrix< double, 6, 6 > m( Matrix< double, 6, 6 >::zeros( ) );
	for ( int i = 0; i < 6; i++ )
		m( i, i ) = 10 - i;

	Vector< double, 6 > s;
	lapack::gesvd( m, s );

	for ( int i = 0; i < 6; i++ )
		BOOST_CHECK_CLOSE( s( i ), m( i, i ), 10e-6 );

	// check symmetric eigenvalue problem
	TestBasicEigenvalues();

	// check Linear Least Squares solver
	TestLinearLeastSquares();

#endif // HAVE_LAPACK
}

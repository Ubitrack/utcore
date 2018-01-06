
#include <utMath/Blas3.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Matrix.h>

#include <typeinfo>
#include <algorithm> //std::generate_n

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <utUtil/BlockTimer.h>
#include <log4cpp/Category.hh>
static log4cpp::Category& timeLogger( log4cpp::Category::getInstance( "Ubitrack.Test.Math.Blas3" ) );

using namespace Ubitrack::Math;

template< typename T >
void testBasicMatrixMatrixProductFunctors( const std::size_t n, const T epsilon )
{	
	
	{
		// random vector parameters
		static const std::size_t size1 = 5;
		static const std::size_t size2 = 3;
		static const std::size_t size3 = size2; // rows should equal columns of first matrix
		static const std::size_t size4 = 4;
		const T max_value( 0.5 );
		typename Random::Matrix< T, size1, size2 >::Uniform randMatrix1( -max_value, max_value );
		typename Random::Matrix< T, size3, size4 >::Uniform randMatrix2( -max_value, max_value );

		// generate some random vectors
		std::vector< Matrix< T, size1, size2 > > matrices1;
		matrices1.reserve( n );
		std::generate_n ( std::back_inserter( matrices1 ), n,  randMatrix1 );
		
		
		std::vector< Matrix< T, size3, size4 > > matrices2;
		matrices2.reserve( n );
		std::generate_n ( std::back_inserter( matrices2 ), n,  randMatrix2 );
		
		
		
		std::stringstream sstream;
		sstream << "blas level 3: " << n << " matrix-matrix multiplications ([" << size1 << "x" << size2 << "]x[" << size3 << "x" << size4 << "]=[" << size1 << "x" << size4 << "])";
		sstream << " of type \"" << typeid( T ).name() << "\" matrices.";
		Ubitrack::Util::BlockTimer g_blas3_mat_mat_product( sstream.str(), timeLogger );
		
		std::vector< Matrix< T, size1, size4 > > results;
		results.reserve( n );

		{
			UBITRACK_TIME( g_blas3_mat_mat_product );
			product( matrices1.begin(), matrices1.end(), matrices2.begin(), std::back_inserter( results ) );
		}
		
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const Matrix< T, size1, size4 > prod = boost::numeric::ublas::prod( matrices1[ i ], matrices2[ i ] );
			T diff = matrixDiff( prod, results[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );			
			//results.push_back( prod );
		}
		
	}
}

void TestBlas3()
{
	testBasicMatrixMatrixProductFunctors< double >( 1000000, 1e-12 );
	testBasicMatrixMatrixProductFunctors< float >( 1000000, 1e-05f );
	
	BOOST_CHECK( 0 == 0 ); // helps to reduce warnings during test phase
}

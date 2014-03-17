
#include <utMath/Blas2.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Matrix.h>


#include <algorithm> //std::generate_n


#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack::Math;


template< typename T >
void testBasicOuterProductFunctors( const std::size_t n, const T epsilon )
{

	{
		// random vector parameters
		static const std::size_t size = 10;
		const T max_value( 0.5 );
		typename Random::Vector< T, size >::Uniform randVector( -max_value, max_value ); 

		// generate some random vectors
		std::vector< Vector< T, size > > vectors;
		vectors.reserve( n );
		std::generate_n ( std::back_inserter( vectors ), n,  randVector );
		
		std::vector< Matrix< T, size, size > > results;
		results.reserve( n );
		outer_product( vectors.begin(), vectors.end(), vectors.begin(), std::back_inserter( results ) );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const Matrix< T, size, size > outer = boost::numeric::ublas::outer_prod( vectors[ i ], vectors[ i ] );
			// results.push_back( outer );
			T diff = matrixDiff( outer, results[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );
		}
	}

	// some 3d and 5d values
	{
		// random vector parameters
		static const std::size_t size1 = 3;
		static const std::size_t size2 = 5;
		const T max_value( 0.5 );
		typename Random::Vector< T, size1 >::Uniform randVector1( -max_value, max_value ); 
		typename Random::Vector< T, size2 >::Uniform randVector2( -max_value, max_value ); 

		// generate some random vectors
		std::vector< Vector< T, size1 > > vectors1;
		vectors1.reserve( n );
		std::generate_n ( std::back_inserter( vectors1 ), n,  randVector1 );
		
		
		// generate some random vectors
		std::vector< Vector< T, size2 > > vectors2;
		vectors2.reserve( n );
		std::generate_n ( std::back_inserter( vectors2 ), n,  randVector2 );
		
		std::vector< Matrix< T, size1, size2 > > results;
		results.reserve( n );
		outer_product( vectors1.begin(), vectors1.end(), vectors2.begin(), std::back_inserter( results ) );

		for( std::size_t i = 0; i<n; ++i )
		{
			const Matrix< T, size1, size2 > outer = boost::numeric::ublas::outer_prod( vectors1[ i ], vectors2[ i ] );
			// results.push_back( outer );
			T diff = matrixDiff( outer, results[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );
		}					 
	}

}


template< typename T >
void testBasiMatrixVectorProductFunctors( const std::size_t n, const T epsilon )
{	
		
	{
		// random vector parameters
		static const std::size_t size1 = 5;
		static const std::size_t size2 = 10;
		const T max_value( 0.5 );
		typename Random::Vector< T, size2 >::Uniform randVector( -max_value, max_value ); 
		typename Random::Matrix< T, size1, size2 >::Uniform randMatrix( -max_value, max_value ); 

		// generate some random vectors
		std::vector< Vector< T, size2 > > vectors;
		vectors.reserve( n );
		std::generate_n ( std::back_inserter( vectors ), n,  randVector );
		
		
		std::vector< Matrix< T, size1, size2 > > matrices;
		matrices.reserve( n );
		std::generate_n ( std::back_inserter( matrices ), n,  randMatrix );
		
		
		std::vector< Vector< T, size1 > > results;
		results.reserve( n );
		product( matrices.begin(), matrices.end(), vectors.begin(), std::back_inserter( results ) );
		
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const Vector< T, size1 > prod = boost::numeric::ublas::prod( matrices[ i ], vectors[ i ] );
			// results.push_back( prod );
			T diff = vectorDiffSum( prod, results[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );
		}
	}
}

void TestBlas2()
{
	
	testBasicOuterProductFunctors< double >( 1000000, 1e-10 );
	testBasiMatrixVectorProductFunctors< double >( 1000000, 1e-12 );
	
	testBasicOuterProductFunctors< float >( 100000, 1e-05f );
	testBasiMatrixVectorProductFunctors< float >( 1000000, 1e-05f );
	
	// BOOST_CHECK( 0 == 0 );
}

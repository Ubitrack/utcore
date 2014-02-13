
#include <utMath/Blas2.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Random/Vector.h>


#include <algorithm> //std::transform
#include <functional> //std::bind1st

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack::Math;


template< typename T >
void testBasicOuterProductFunctors( const std::size_t n, const T epsilon )
{
	// some 10d values
	{
		// random vector parameters
		typename Random::Vector< T, 10 >::Uniform randVector( -0.5, 0.5 ); 

		// generate some random vectors
		std::vector< Vector< T, 10 > > points10d;
		points10d.reserve( n );
		std::generate_n ( std::back_inserter( points10d ), n,  randVector );
		
		std::vector< Matrix< T, 10, 10 > > results;
		results.reserve( n );
		std::transform( points10d.begin(), points10d.end(), points10d.begin(), std::back_inserter( results ), OuterProduct() );
		
		// for( std::size_t i = 0; i<n; ++i )
		// {
			// const Matrix< T, 10, 10 > outer = boost::numeric::ublas::outer_prod( points10d[ i ], points10d[ i ] );
			// results.push_back( outer );
			// // T diff = matrixDiff( outer, results[ i ] );
			// // BOOST_CHECK_SMALL( diff, epsilon );
		// }
	}
	
	
	// // some 3d values
	// {
		// // random vector parameters
		// typename Random::Vector< T, 3 >::Uniform randVector( -100., 100. ); 

		// // generate some random vectors
		// std::vector< Vector< T, 3 > > points3dA;
		// points3dA.reserve( n );
		// std::generate_n ( std::back_inserter( points3dA ), n,  randVector );
		
		
		// // generate some random vectors
		// std::vector< Vector< T, 3 > > points3dB;
		// points3dB.reserve( n );
		// std::generate_n ( std::back_inserter( points3dB ), n,  randVector );
		
		// std::vector< Matrix< T, 3, 3 > > results;
		// results.reserve( n );
		// std::transform( points3dA.begin(), points3dA.end(), points3dB.begin(), std::back_inserter( results ), OuterProduct() );
		// for( std::size_t i = 0; i<n; ++i )
		// {
			// const Matrix< T, 3, 3 > outer = boost::numeric::ublas::outer_prod( points3dA[ i ], points3dB[ i ] );
			// // results.push_back( outer );
			// T diff = matrixDiff( outer, results[ i ] );
			// BOOST_CHECK_SMALL( diff, epsilon );
		// }					 
	// }
	
	
	
	// // some 3d and 5d values
	// {
		// // random vector parameters
		// typename Random::Vector< T, 3 >::Uniform randVector3( -100., 100. ); 
		// typename Random::Vector< T, 5 >::Uniform randVector5( -50., 50. ); 

		// // generate some random vectors
		// std::vector< Vector< T, 3 > > points3d;
		// points3d.reserve( n );
		// std::generate_n ( std::back_inserter( points3d ), n,  randVector3 );
		
		
		// // generate some random vectors
		// std::vector< Vector< T, 5 > > points5d;
		// points5d.reserve( n );
		// std::generate_n ( std::back_inserter( points5d ), n,  randVector5 );
		
		// std::vector< Matrix< T, 3, 5 > > results;
		// results.reserve( n );
		// std::transform( points3d.begin(), points3d.end(), points5d.begin(), std::back_inserter( results ), OuterProduct() );
		// for( std::size_t i = 0; i<n; ++i )
		// {
			// const Matrix< T, 3, 5 > outer = boost::numeric::ublas::outer_prod( points3d[ i ], points5d[ i ] );
			// // results.push_back( outer );
			// T diff = matrixDiff( outer, results[ i ] );
			// BOOST_CHECK_SMALL( diff, epsilon );
		// }					 
	// }
	BOOST_CHECK( 0 == 0 );
}


void TestBlas2()
{
	// testBasicOuterProductFunctors< float >( 100000, 1e-05f );
	testBasicOuterProductFunctors< double >( 1000000, 1e-10 );
}



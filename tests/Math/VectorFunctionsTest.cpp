
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Blas1.h>
#include <utMath/VectorFunctions.h>
#include <utMath/Random/Vector.h>

#include <algorithm> //std::transform

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack::Math;



template< typename T >
void testDistance( const std::size_t n, const T epsilon )
{
	// some 10d values
	{
		typedef Vector< T, 10 > vector_type;
		
		// random vector parameters
		typename Random::Vector< T, 10 >::Uniform randVector( -50., 50. ); 

		// generate some random vectors
		std::vector< vector_type > points10d;
		points10d.reserve( n );
		std::generate_n ( std::back_inserter( points10d ), n, randVector );
		
		std::vector< T > results;
		results.reserve( n );
		std::transform( points10d.begin(), points10d.end(), std::back_inserter( results ), Distance< vector_type >() );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T dist = results[ i ] - norm_2( points10d[ i ] );
			BOOST_CHECK_SMALL( dist, epsilon );
		}
	}
	
	// some 3d values
	{
	
		typedef Vector< T, 3 > vector_type;
		// random vector parameters
		typename Random::Vector< T, 3 >::Uniform randVector( -100., 100. ); 

		// generate some random vectors
		std::vector< vector_type > points3dA;
		points3dA.reserve( n );
		std::generate_n ( std::back_inserter( points3dA ), n,  randVector );
		
		
		// generate some random vectors
		std::vector< vector_type > points3dB;
		points3dB.reserve( n );
		std::generate_n ( std::back_inserter( points3dB ), n, randVector );
		
		std::vector< T > results;
		results.reserve( n );
		std::transform( points3dA.begin(), points3dA.end(), points3dB.begin(), std::back_inserter( results ), Distance< vector_type >() );
		for( std::size_t i = 0; i<n; ++i )
		{
			const vector_type diffVector = points3dA[ i ] - points3dB[ i ];
			const T dist = results[ i ] - norm_2( diffVector );
			BOOST_CHECK_SMALL( dist, epsilon );
		}					 
	}
	
	// some 3d values
	{
	
		typedef Vector< T, 3 > vector_type;
		// random vector parameters
		typename Random::Vector< T, 3 >::Uniform randVector( -100., 100. ); 

		// generate some random vectors
		std::vector< vector_type > points3dA;
		points3dA.reserve( n );
		std::generate_n ( std::back_inserter( points3dA ), n,  randVector );
		
		
		// generate some random vectors
		std::vector< vector_type > points3dB;
		points3dB.reserve( n );
		std::generate_n ( std::back_inserter( points3dB ), n, randVector );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T result = distance( points3dA[ i ], points3dB[ i ] );
			const vector_type diffVector = points3dA[ i ] - points3dB[ i ];
			const T dist = result - norm_2( diffVector );
			BOOST_CHECK_SMALL( dist, epsilon );
		}					 
	}
}

template< typename T >
void testNormalization( const std::size_t n, const T epsilon )
{
	// some 10d values
	{
		typedef Vector< T, 10 > vector_type;
		// random vector parameters
		typename Random::Vector< T, 10 >::Uniform randVector( -0.5, 0.5 ); 

		// generate some random vectors
		std::vector< vector_type > points10d;
		points10d.reserve( n );
		std::generate_n ( std::back_inserter( points10d ), n, randVector );
		
		std::vector< vector_type > results;
		results.reserve( n );
		std::transform( points10d.begin(), points10d.end(), std::back_inserter( results ), Normalize< vector_type >() );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T diff = 1. - norm_2( results[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );
		}
	}
	
	
	// some 3d values
	{
		typedef Vector< T, 3 > vector_type;
		// random vector parameters
		typename Random::Vector< T, 3 >::Uniform randVector( -100., 100. ); 

		// generate some random vectors
		std::vector< vector_type > points3dA;
		points3dA.reserve( n );
		std::generate_n ( std::back_inserter( points3dA ), n,  randVector );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const vector_type result = normalize( points3dA[ i ] );
			const T diff = 1. - norm_2( result );
			BOOST_CHECK_SMALL( diff, epsilon );
		}					 
	}
}


void TestVectorFunctions()
{
	testDistance< double >( 1000000, 1e-10 );
	testNormalization< double >( 1000000, 1e-10 );
	
	testDistance< float >( 1000000, 1e-10f );
	testNormalization< float >( 1000000, 1e-06f );
}




#include <utMath/Blas1.h>
#include <utMath/Vector.h>
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>


#include <algorithm> //std::transform
#include <functional> //std::bind1st

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack::Math;


template< typename T >
void testBasicInnerProductFunctors( const std::size_t n, const T epsilon )
{
	// some 10d values
	{
		// random vector parameters
		typename Random::Vector< T, 10 >::Uniform randVector( -0.5, 0.5 ); 

		// generate some random vectors
		std::vector< Vector< T, 10 > > points10d;
		points10d.reserve( n );
		std::generate_n ( std::back_inserter( points10d ), n,  randVector );
		
		std::vector< T > results;
		results.reserve( n );
		std::transform( points10d.begin(), points10d.end(), std::back_inserter( results ), InnerProduct() );
		for( std::size_t i = 0; i<n; ++i )
		{
			const T inner = boost::numeric::ublas::inner_prod( points10d[ i ], points10d[ i ] );
			// results.push_back( inner );
			BOOST_CHECK_SMALL( inner - results[ i ], epsilon );//inner << " vs. " << results[ i ] );
		}
	}
	
	
	// some 3d values
	{
		// random vector parameters
		typename Random::Vector< T, 3 >::Uniform randVector( -100., 100. ); 

		// generate some random vectors
		std::vector< Vector< T, 3 > > points3dA;
		points3dA.reserve( n );
		std::generate_n ( std::back_inserter( points3dA ), n,  randVector );
		
		
		// generate some random vectors
		std::vector< Vector< T, 3 > > points3dB;
		points3dB.reserve( n );
		std::generate_n ( std::back_inserter( points3dB ), n,  randVector );
		
		std::vector< T > results;
		results.reserve( n );
		std::transform( points3dA.begin(), points3dA.end(), points3dB.begin(), std::back_inserter( results ), InnerProduct() );
		for( std::size_t i = 0; i<n; ++i )
		{
			const T inner = boost::numeric::ublas::inner_prod( points3dA[ i ], points3dB[ i ] );
			// results.push_back( inner );
			BOOST_CHECK_SMALL( inner - results[ i ], epsilon );
		}
	}
}


template< typename T >
void testBasicNorm2Functors( const std::size_t n, const T epsilon )
{
	// some 10d values
	{
		// random vector parameters
		typename Random::Vector< T, 10 >::Uniform randVector( -20., 20. ); 

		// generate some random vectors
		std::vector< Vector< T, 10 > > points10d;
		points10d.reserve( n );
		std::generate_n ( std::back_inserter( points10d ), n,  randVector );
		
		std::vector< T > results;
		results.reserve( n );
		std::transform( points10d.begin(), points10d.end(), std::back_inserter( results ), Norm_2() );
		for( std::size_t i = 0; i<n; ++i )
		{
			const T inner = boost::numeric::ublas::norm_2( points10d[ i ] );
			// results.push_back( inner );
			BOOST_CHECK_SMALL( inner - results[ i ], epsilon );//inner << " vs. " << results[ i ] );
		}
	}
	
	
	// some 3d values
	{
		// random vector parameters
		typename Random::Vector< T, 3 >::Uniform randVector( -100., 100. ); 

		// generate some random vectors
		std::vector< Vector< T, 3 > > points3d;
		points3d.reserve( n );
		std::generate_n ( std::back_inserter( points3d ), n,  randVector );
		
		
		std::vector< T > results;
		results.reserve( n );
		std::transform( points3d.begin(), points3d.end(), std::back_inserter( results ), Norm_2() );
		for( std::size_t i = 0; i<n; ++i )
		{
			const T inner = boost::numeric::ublas::norm_2( points3d[ i ] );
			// results.push_back( inner );
			BOOST_CHECK_SMALL( inner - results[ i ], epsilon );
		}
	}
}


template< typename T >
void testBasicNorm1Functors( const std::size_t n, const T epsilon )
{
	// some 10d values
	{
		// random vector parameters
		typename Random::Vector< T, 10 >::Uniform randVector( -20., 20. ); 

		// generate some random vectors
		std::vector< Vector< T, 10 > > points10d;
		points10d.reserve( n );
		std::generate_n ( std::back_inserter( points10d ), n,  randVector );
		
		std::vector< T > results;
		results.reserve( n );
		std::transform( points10d.begin(), points10d.end(), std::back_inserter( results ), Norm_1() );
		for( std::size_t i = 0; i<n; ++i )
		{
			const T inner = boost::numeric::ublas::norm_1( points10d[ i ] );
			// results.push_back( inner );
			BOOST_CHECK_SMALL( inner - results[ i ], epsilon );//inner << " vs. " << results[ i ] );
		}
	}
}

void TestBlas1()
{
	testBasicInnerProductFunctors< float >( 100000, 1e-02f );
	testBasicInnerProductFunctors< double >( 100000, 1e-10 );
	testBasicNorm2Functors< float >( 100000, 1e-02f );
	testBasicNorm2Functors< double >( 1000000, 1e-10 );
	testBasicNorm1Functors< float >( 100000, 1e-02f );
	testBasicNorm1Functors< double >( 1000000, 1e-10 );
}



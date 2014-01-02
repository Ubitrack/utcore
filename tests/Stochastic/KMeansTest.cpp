
#include <utMath/Blas1.h> //-> norm_2
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Stochastic/k_means.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack::Math;

template< typename T >
void testBasicKMeans( const std::size_t n_runs, const std::size_t max_n, const std::size_t max_cluster )
{

	for( std::size_t run = 0; run< n_runs; ++run )
	{
		typename Random::Vector< T, 2 >::Uniform randPoints2D( -5, 5 ); 
		
		const std::size_t cluster = Random::distribute_uniform( static_cast< std::size_t > ( 3 ), max_cluster );
		const std::size_t n = Random::distribute_uniform( max_n/2, max_n );
		
		// generate some random points
		std::vector< Vector< T, 2 > > pts2D;
		pts2D.reserve( n );
		std::generate_n ( std::back_inserter( pts2D ), n,  randPoints2D );
		
		//reserve some space for indices
		std::vector< std::size_t > indices;
		indices.reserve( n );
		
		//reserve some space for mean values 
		std::vector< Vector< T, 2 > > centroids;
		centroids.reserve( n );
		
		// execute the k-means algorithm
		Stochastic::k_means( pts2D.begin(), pts2D.end(), cluster, std::back_inserter( centroids ), std::back_inserter( indices ) );
		
		// std::cout << "Centroids:\n " << centroids << "\n";
		
		for( std::size_t i = 0; i<n; ++i)
		{
			const std::size_t index = indices[ i ];
			Vector< T, 2 > diffVec0 = pts2D[ i ] - centroids[ index ];
			const T diff0 = Norm_2()( diffVec0 );
			
			for( std::size_t k = 0; k<cluster; ++k)
			{
				if( index == k )
					continue;
				Vector< T, 2 > diffVec1 = pts2D[ i ] - centroids[ k ];
				const T diff1 = Norm_2()( diffVec1 );
				BOOST_CHECK_MESSAGE( diff0 <= diff1, cluster << " cluster and " << n << " values, remaining difference: " << diff0 << " vs. " << diff1 << "." );
			}
		}
	}
}

void TestKMeans()
{
	testBasicKMeans< double >( 10, 10000, 5 );
	testBasicKMeans< float >( 10, 10000, 5 );
}


#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Geometry/Conic.h>

#include <algorithm>

#include "../tools.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack::Math;


template< typename T >
void testBasicConicFunctors( const std::size_t n )
{
	// random conic parameters
	typename Random::Vector< 6, T >::Uniform randConic( -0.5, 0.5 ); 

	// generate some random conics
	std::vector< Ubitrack::Math::Vector< 6, T > > conics;
	conics.reserve( n );
	std::generate_n ( std::back_inserter( conics ), n,  randConic );

	// generate matrix representation
	std::vector< Ubitrack::Math::Matrix< 3, 3, T > > conic_matrices;
	conic_matrices.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_matrices ), Geometry::matrix_from_conic< T >() );

	// generate vector representation, again
	std::vector< Ubitrack::Math::Vector< 6, T > > conics_again;
	conics_again.reserve( n );
	std::transform( conic_matrices.begin(), conic_matrices.end(), std::back_inserter( conics_again ), Geometry::conic_from_matrix< T >() );

	//check if conic is still the same
	for( std::size_t i( 0 ); i<n; ++i )
		BOOST_CHECK_SMALL( vectorDiff( conics_again[ i ], conics[ i ] ), 10e-20 );

	// invert conics
	std::vector< Ubitrack::Math::Vector< 6, T > > inv_conics;
	inv_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( inv_conics ), Geometry::conic_inverse< T >() );
	
	// calculate conic determinants
	std::vector< T > conic_determinant;
	conic_determinant.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_determinant ), Geometry::conic_determinant< T >() );

	// estimate conic semi-axes
	std::vector< Ubitrack::Math::Vector< 2, T > > semi_axes;
	semi_axes.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( semi_axes ), Geometry::conic_semi_axes< T >() );

	// estimate conic angles
	std::vector< T > angles;
	angles.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( angles ), Geometry::conic_angle< T >() );

	// estimate conic center
	std::vector< Ubitrack::Math::Vector< 2, T > > centers;
	centers.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( centers ), Geometry::conic_center< T >() );

	// estimate conic eccentricity
	std::vector< T > eccentricities;
	eccentricities.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( eccentricities ), Geometry::conic_eccentricity< T >() );

	// calculate conic determinantes
	std::vector< T > conic_areas;
	conic_areas.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_areas ), Geometry::conic_area< T >() );

	// scale the conics
	std::vector< Ubitrack::Math::Vector< 6, T > > scaled_conics;
	scaled_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( scaled_conics ), Geometry::scale_conic_unsafe< T >( 10.0 ) );

	// translate the conics
	/*std::vector< Ubitrack::Math::Vector< 6, T > > translated_conics;
	translated_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( translated_conics ), Geometry::translate_conic< T >() );
	*/

	// next steps are even more useless, just check if they compile
	std::size_t n_c = std::count_if( conics.begin(), conics.end(), Geometry::is_conic_circle< T >() );
	std::size_t n_d = std::count_if( conics.begin(), conics.end(), Geometry::is_conic_degenerate< T >() );
	std::size_t n_e = std::count_if( conics.begin(), conics.end(), Geometry::is_conic_ellipse< T >() );
	std::size_t n_p = std::count_if( conics.begin(), conics.end(), Geometry::is_conic_parabola< T >() );

	//std::cout << "Conics " << n_c << " " << n_d << " " << n_e << " " << n_p << std::endl;

	// estimate conics' upper and lower limit
	std::vector< Ubitrack::Math::Vector< 2, T > > conic_ull;
	conic_ull.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_ull ), Geometry::conic_upper_lower_limit< T >() );

	// estimate conics' left and right limit
	std::vector< Ubitrack::Math::Vector< 2, T > > conic_lrl;
	conic_lrl.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_lrl ), Geometry::conic_left_right_limit< T >() );

}


void TestConic()
{
	testBasicConicFunctors< float >( 10000 );
	testBasicConicFunctors< double >( 10000 );
}



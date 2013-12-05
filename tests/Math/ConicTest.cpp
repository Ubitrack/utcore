
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include <utMath/Geometry/Conic.h>
#include <utMath/Geometry/QuadricFunctors.h>
#include <utMath/Functors/VectorFunctors.h>

#include <algorithm> //std::transform
#include <functional> //std::bind1st

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
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_matrices ), Geometry::MatrixFromConic< T >() );

	// generate vector representation, again
	std::vector< Ubitrack::Math::Vector< 6, T > > conics_again;
	conics_again.reserve( n );
	std::transform( conic_matrices.begin(), conic_matrices.end(), std::back_inserter( conics_again ), Geometry::ConicFromMatrix< T >() );

	//check if conic is still the same
	for( std::size_t i( 0 ); i<n; ++i )
		BOOST_CHECK_SMALL( vectorDiff( conics_again[ i ], conics[ i ] ), static_cast< T > ( 10e-20 ) );

	// invert conics
	std::vector< Ubitrack::Math::Vector< 6, T > > inv_conics;
	inv_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( inv_conics ), Geometry::ConicInverse< T >() );
	
	// calculate conic determinants
	std::vector< T > conic_determinant;
	conic_determinant.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_determinant ), Geometry::ConicDeterminant< T >() );

	// estimate conic semi-axes
	std::vector< Ubitrack::Math::Vector< 2, T > > semi_axes;
	semi_axes.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( semi_axes ), Geometry::ConicSemiAxes< T >() );

	// estimate conic angles
	std::vector< T > angles;
	angles.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( angles ), Geometry::ConicAngle< T >() );

	// estimate conic center
	std::vector< Ubitrack::Math::Vector< 2, T > > centers;
	centers.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( centers ), Geometry::ConicCenter< T >() );

	// estimate conic eccentricity
	std::vector< T > eccentricities;
	eccentricities.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( eccentricities ), Geometry::ConicEccentricity< T >() );

	// calculate conic determinantes
	std::vector< T > conic_areas;
	conic_areas.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_areas ), Geometry::ConicArea< T >() );

	// scale the conics
	std::vector< Ubitrack::Math::Vector< 6, T > > scaled_conics;
	scaled_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( scaled_conics ), Geometry::ScaleConicUnsafe< T >( 10.0 ) );

	// translate the conics
	/*std::vector< Ubitrack::Math::Vector< 6, T > > translated_conics;
	translated_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( translated_conics ), Geometry::TranslateConic< T >() );
	*/

	// next steps are even more useless, just check if they compile
	std::size_t n_c = std::count_if( conics.begin(), conics.end(), Geometry::IsConicCircle< T >() );
	std::size_t n_d = std::count_if( conics.begin(), conics.end(), Geometry::IsConicDegenerate< T >() );
	std::size_t n_e = std::count_if( conics.begin(), conics.end(), Geometry::IsConicEllipse< T >() );
	std::size_t n_p = std::count_if( conics.begin(), conics.end(), Geometry::IsConicParabola< T >() );

	std::cout << "Conics " << n_c << " " << n_d << " " << n_e << " " << n_p << std::endl;

	// estimate conics' upper and lower limit
	std::vector< Ubitrack::Math::Vector< 2, T > > conic_ull;
	conic_ull.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_ull ), Geometry::ConicUpperLowerLimit< T >() );

	// estimate conics' left and right limit
	std::vector< Ubitrack::Math::Vector< 2, T > > conic_lrl;
	conic_lrl.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_lrl ), Geometry::ConicLeftRightLimit< T >() );

}



template< typename T >
void testRandomQuadricProjection( const std::size_t n )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< 3, T >::Uniform randTranslation( -100, 100 );
	Pose camPose( randQuat() , randTranslation() );
	Matrix< 3, 4, T > projection( camPose );
	
	typename Random::Vector< 4, T >::Uniform randSpheroid( -5.0, 5.0 );
	typename Random::Vector< 6, T >::Uniform randEllipsoid( -5.0, 5.0 );

	// generate some random quadrics
	std::vector< Ubitrack::Math::Vector< 6, T > > ellipsoids;
	ellipsoids.reserve( n );
	std::generate_n ( std::back_inserter( ellipsoids ), n,  randEllipsoid );

	//std::copy( ellipsoids.begin(), ellipsoids.end(), std::ostream_iterator< Ubitrack::Math::Vector< 6, T > > ( std::cout, ", ") );

	// project these ellipsoids
	std::vector< Ubitrack::Math::Vector< 6, T > > conics1;
	conics1.reserve( n );
	std::transform( ellipsoids.begin(), ellipsoids.end(), std::back_inserter( conics1 ), std::bind1st( Geometry::ProjectEllipsoid< T >(), projection ) );
	//std::transform( conics1.begin(), conics1.end(), conics1.begin(), Ubitrack::Math::Functors::NormalizeVector< 6, T >() );

	//generate some quadrics
	std::vector< Ubitrack::Math::Vector< 10, T > > quadrics1;
	quadrics1.reserve( n );
	std::transform( ellipsoids.begin(), ellipsoids.end(), std::back_inserter( quadrics1 ), Geometry::Ellipsoid2Quadric< T >() );

	// project these quadrics
	std::vector< Ubitrack::Math::Vector< 6, T > > conics2;
	conics2.reserve( n );
	std::transform( quadrics1.begin(), quadrics1.end(), std::back_inserter( conics2 ), std::bind1st( Geometry::ProjectQuadric< T >(), projection ) );
	//std::transform( conics2.begin(), conics2.end(), conics2.begin(), Ubitrack::Math::Functors::NormalizeVector< 6, T >() );


	// write out conics from ellipsoids and quadrics
	/*std::cout << std::endl;
	std::copy( conics1.begin(), conics1.end(), std::ostream_iterator< Ubitrack::Math::Vector< 6, T > > ( std::cout, ", ") );
	std::cout << std::endl;
	std::copy( conics2.begin(), conics2.end(), std::ostream_iterator< Ubitrack::Math::Vector< 6, T > > ( std::cout, ", ") );*/

	//another test

	// generate some random spheroids
	std::vector< Ubitrack::Math::Vector< 4, T > > spheroids;
	spheroids.reserve( n );
	std::generate_n ( std::back_inserter( spheroids ), n,  randSpheroid );

	
	// project the spheroids
	std::vector< Ubitrack::Math::Vector< 6, T > > conics3;
	conics3.reserve( n );
	std::transform( spheroids.begin(), spheroids.end(), std::back_inserter( conics3 ), std::bind1st( Geometry::ProjectSpheroid< T >(), projection ) );

	/*
	// commented out, since it does not really represent the same values.
	// generate some corresponding quadrics
	std::vector< Ubitrack::Math::Vector< 10, T > > quadrics2;

	for( std::size_t i( 0 ); i<n; ++i )
	{
		Ubitrack::Math::Vector< 10, T > quadric;

		quadric( 0 ) = spheroids[ i ]( 0 );
		quadric( 1 ) = spheroids[ i ]( 0 );
		quadric( 2 ) = spheroids[ i ]( 1 );
		quadric( 3 ) = quadric( 4 ) = quadric( 5 ) = 0;
		quadric( 6 ) = spheroids[ i ]( 2 );
		quadric( 7 ) = spheroids[ i ]( 3 );
		quadric( 8 ) = spheroids[ i ]( 1 );
		quadric( 9 ) = -1;
		quadrics2.push_back( quadric );
	}

	// project the quadrics
	std::vector< Ubitrack::Math::Vector< 6, T > > conics4;
	conics4.reserve( n );
	std::transform( quadrics2.begin(), quadrics2.end(), std::back_inserter( conics4 ), std::bind1st( Geometry::ProjectQuadric< T >(), projection ) );
	*/

	// write out conics from spheroids and quadrics
	/*std::cout << std::endl;
	std::copy( conics3.begin(), conics3.end(), std::ostream_iterator< Ubitrack::Math::Vector< 6, T > > ( std::cout, ", ") );
	std::cout << std::endl;
	std::copy( conics4.begin(), conics4.end(), std::ostream_iterator< Ubitrack::Math::Vector< 6, T > > ( std::cout, ", ") );*/

}


void TestConic()
{
	// float is usually not sufficient here
	testBasicConicFunctors< float >( 10000 );
	testBasicConicFunctors< double >( 10000 );
	testRandomQuadricProjection< float >( 10000 );
	testRandomQuadricProjection< double >( 10000 );
}



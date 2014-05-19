
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/VectorFunctions.h>
#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>
#include <utMath/Geometry/Conic.h>
#include <utMath/Geometry/ConicCovariance.h>

#include <utMath/Geometry/QuadricFunctors.h>


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
	typename Random::Vector< T, 6 >::Uniform randConic( -0.5, 0.5 ); 

	// generate some random conics
	std::vector< Ubitrack::Math::Vector< T, 6 > > conics;
	conics.reserve( n );
	std::generate_n ( std::back_inserter( conics ), n,  randConic );

	// generate matrix representation
	std::vector< Ubitrack::Math::Matrix< T, 3, 3 > > conic_matrices;
	conic_matrices.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_matrices ), Geometry::MatrixFromConic< T >() );

	// generate vector representation, again
	std::vector< Ubitrack::Math::Vector< T, 6 > > conics_again;
	conics_again.reserve( n );
	std::transform( conic_matrices.begin(), conic_matrices.end(), std::back_inserter( conics_again ), Geometry::ConicFromMatrix< T >() );

	//check if conic is still the same
	for( std::size_t i( 0 ); i<n; ++i )
		BOOST_CHECK_SMALL( vectorDiff( conics_again[ i ], conics[ i ] ), static_cast< T > ( 10e-20 ) );

	// invert conics
	std::vector< Ubitrack::Math::Vector< T, 6 > > inv_conics;
	inv_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( inv_conics ), Geometry::ConicInverse< T >() );
	
	// calculate conic determinants
	std::vector< T > conic_determinant;
	conic_determinant.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_determinant ), Geometry::ConicDeterminant< T >() );

	// estimate conic semi-axes
	std::vector< Ubitrack::Math::Vector< T, 2 > > semi_axes;
	semi_axes.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( semi_axes ), Geometry::ConicSemiAxes< T >() );

	// estimate conic angles
	std::vector< T > angles;
	angles.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( angles ), Geometry::ConicAngle< T >() );

	// estimate conic center
	std::vector< Ubitrack::Math::Vector< T, 2 > > centers;
	centers.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( centers ), Geometry::ConicCenter< T >() );

	// estimate conic eccentricity
	std::vector< T > eccentricities;
	eccentricities.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( eccentricities ), Geometry::ConicEccentricity< T >() );

	// calculate conic determinants
	std::vector< T > conic_areas;
	conic_areas.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_areas ), Geometry::ConicArea< T >() );

	// scale the conics
	std::vector< Ubitrack::Math::Vector< T, 6 > > scaled_conics;
	scaled_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( scaled_conics ), Geometry::ScaleConicUnsafe< T >( 10.0 ) );

	// translate the conics
	/*std::vector< Ubitrack::Math::Vector< T, 6 > > translated_conics;
	translated_conics.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( translated_conics ), Geometry::TranslateConic< T >() );
	*/

	
	// estimate conics' upper and lower limit
	std::vector< Ubitrack::Math::Vector< T, 2 > > conic_ull;
	conic_ull.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_ull ), Geometry::ConicUpperLowerLimit< T >() );

	// estimate conics' left and right limit
	std::vector< Ubitrack::Math::Vector< T, 2 > > conic_lrl;
	conic_lrl.reserve( n );
	std::transform( conics.begin(), conics.end(), std::back_inserter( conic_lrl ), Geometry::ConicLeftRightLimit< T >() );
	
	
	{	// covariance estimation of point conic
		std::vector< Ubitrack::Math::Vector< T, 2 > >::const_iterator it = conic_ull.begin();
		std::vector< Ubitrack::Math::Vector< T, 6 > >::const_iterator itC = conics.begin();
		
		for( ; it != conic_ull.end(); ++it, ++itC )
		{
			const T dist = std::fabs( (*it)[ 0 ] - (*it)[ 1 ]);
			int start = std::ceil( (*it)[ 0 ] );
			int end = std::floor( (*it)[ 1 ] );
			if( dist != dist || dist < 5 || !Geometry::IsConicEllipse< T >()( *itC ) ) // only apply to ellipses
				continue;
			
			// std::cout << " Elements " << std::fabs( (*it)[ 0 ] - (*it)[ 1 ]) << std::endl;			
			
			std::vector< Ubitrack::Math::Vector< T, 2 > > pixels;
			// pixels.reserve( );
			
			for( int y = start; y < end; ++y )
			{
				Ubitrack::Math::Vector< T, 2 > x12 = Geometry::ConicHorizontalIntersection< T >()( *itC, y );
				// std::cout << "Vec y " << y << std::endl;
				pixels.push_back( Ubitrack::Math::Vector< T, 2 >( x12[ 0 ], y ) );
				pixels.push_back( Ubitrack::Math::Vector< T, 2 >( x12[ 1 ], y ) );
			}
			Ubitrack::Math::Matrix< T, 6, 6 > cov;
			Ubitrack::Math::Geometry::estimateCovariance(  pixels.begin(), pixels.end(), *itC, cov );
			// std::cout << "Conic      " << *itC << std::endl;
			// std::cout << "Covariance " << cov << "\n";
		}
	}
	
	// some conic normalization
	std::transform( conics.begin(), conics.end(), conics.begin(), Ubitrack::Math::Normalize< Ubitrack::Math::Vector< T, 6 > >() );
	
	// next steps are even more useless, just check if they compile
	std::size_t n_c = std::count_if( conics.begin(), conics.end(), Geometry::IsConicCircle< T >() );
	std::size_t n_d = std::count_if( conics.begin(), conics.end(), Geometry::IsConicDegenerate< T >() );
	std::size_t n_e = std::count_if( conics.begin(), conics.end(), Geometry::IsConicEllipse< T >() );
	std::size_t n_h = std::count_if( conics.begin(), conics.end(), Geometry::IsConicHyperbola< T >() );
	std::size_t n_p = std::count_if( conics.begin(), conics.end(), Geometry::IsConicParabola< T >() );

	// std::cout << "Conics " << n_c << " " << n_d << " " << n_e << " " << n_h << " " << n_p << std::endl;
	const std::size_t n_sum = n_c + n_d + n_e + n_h + n_p;
	BOOST_CHECK( n_sum >= n );

}



template< typename T >
void testRandomQuadricProjection( const std::size_t n )
{
	typename Random::Quaternion< T >::Uniform randQuat;
	typename Random::Vector< T, 3 >::Uniform randTranslation( -100, 100 );
	Pose camPose( randQuat() , randTranslation() );
	Matrix< T, 3, 4 > projection( camPose );
	
	typename Random::Vector< T, 4 >::Uniform randSpheroid( -5.0, 5.0 );
	typename Random::Vector< T, 6 >::Uniform randEllipsoid( -5.0, 5.0 );

	// generate some random quadrics
	std::vector< Ubitrack::Math::Vector< T, 6 > > ellipsoids;
	ellipsoids.reserve( n );
	std::generate_n ( std::back_inserter( ellipsoids ), n,  randEllipsoid );

	//std::copy( ellipsoids.begin(), ellipsoids.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 6 > > ( std::cout, ", ") );

	// project these ellipsoids
	std::vector< Ubitrack::Math::Vector< T, 6 > > conics1;
	conics1.reserve( n );
	std::transform( ellipsoids.begin(), ellipsoids.end(), std::back_inserter( conics1 ), std::bind1st( Geometry::ProjectEllipsoid< T >(), projection ) );
	//std::transform( conics1.begin(), conics1.end(), conics1.begin(), Ubitrack::Math::NormalizeVector() );

	//generate some quadrics
	std::vector< Ubitrack::Math::Vector< T, 10 > > quadrics1;
	quadrics1.reserve( n );
	std::transform( ellipsoids.begin(), ellipsoids.end(), std::back_inserter( quadrics1 ), Geometry::Ellipsoid2Quadric< T >() );

	// project these quadrics
	std::vector< Ubitrack::Math::Vector< T, 6 > > conics2;
	conics2.reserve( n );
	std::transform( quadrics1.begin(), quadrics1.end(), std::back_inserter( conics2 ), std::bind1st( Geometry::ProjectQuadric< T >(), projection ) );
	//std::transform( conics2.begin(), conics2.end(), conics2.begin(), Ubitrack::Math::NormalizeVector() );


	// write out conics from ellipsoids and quadrics
	/*std::cout << std::endl;
	std::copy( conics1.begin(), conics1.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 6 > > ( std::cout, ", ") );
	std::cout << std::endl;
	std::copy( conics2.begin(), conics2.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 6 > > ( std::cout, ", ") );*/

	//another test

	// generate some random spheroids
	std::vector< Ubitrack::Math::Vector< T, 4 > > spheroids;
	spheroids.reserve( n );
	std::generate_n ( std::back_inserter( spheroids ), n,  randSpheroid );

	
	// project the spheroids
	std::vector< Ubitrack::Math::Vector< T, 6 > > conics3;
	conics3.reserve( n );
	std::transform( spheroids.begin(), spheroids.end(), std::back_inserter( conics3 ), std::bind1st( Geometry::ProjectSpheroid< T >(), projection ) );

	/*
	// commented out, since it does not really represent the same values.
	// generate some corresponding quadrics
	std::vector< Ubitrack::Math::Vector< T, 10 > > quadrics2;

	for( std::size_t i( 0 ); i<n; ++i )
	{
		Ubitrack::Math::Vector< T, 10 > quadric;

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
	std::vector< Ubitrack::Math::Vector< T, 6 > > conics4;
	conics4.reserve( n );
	std::transform( quadrics2.begin(), quadrics2.end(), std::back_inserter( conics4 ), std::bind1st( Geometry::ProjectQuadric< T >(), projection ) );
	*/

	// write out conics from spheroids and quadrics
	/*std::cout << std::endl;
	std::copy( conics3.begin(), conics3.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 6 > > ( std::cout, ", ") );
	std::cout << std::endl;
	std::copy( conics4.begin(), conics4.end(), std::ostream_iterator< Ubitrack::Math::Vector< T, 6 > > ( std::cout, ", ") );*/

}


void TestConic()
{
	// float is usually not sufficient here
	testBasicConicFunctors< float >( 10000 );
	testBasicConicFunctors< double >( 10000 );
	testRandomQuadricProjection< float >( 10000 );
	testRandomQuadricProjection< double >( 10000 );
}



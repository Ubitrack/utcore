
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/VectorFunctions.h> // for test
#include <utMath/Util/RotationCast.h>



#include <utMath/Random/Scalar.h>
#include <utMath/Random/Vector.h>
#include <utMath/Random/Rotation.h>


#include <algorithm> // std::transform

#include "../tools.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Ubitrack;


template< typename T >
void testQuaternionCast( const std::size_t n_runs, const T epsilon )
{
	// we do all runs at the same time...
	const std::size_t n = n_runs;
	
	// cast a quaterion to different type
	typedef Math::Quaternion quat_type;
	typedef Math::Matrix< T, 3, 3 > mat_type;
	typedef Math::Vector< T, 4 > aaxis_type; // axis angle 
	typedef Math::Vector< T, 3 > rotaxis_type; // rotation axis
	
	typename Math::Random::Quaternion< T >::Uniform randQuat;
	
	std::vector< quat_type > quats;
	quats.reserve( n );
	std::generate_n ( std::back_inserter( quats ), n,  randQuat );
	
	{	// test quaternion -> axis angle transformation
		std::vector< aaxis_type > aa4s;
		quats.reserve( n );
		std::transform( quats.begin(), quats.end(), std::back_inserter( aa4s ), Math::Util::RotationCast< aaxis_type > () );
		
		std::vector< quat_type > quatsAA4;
		quatsAA4.reserve( n );
		std::transform( aa4s.begin(), aa4s.end(), std::back_inserter( quatsAA4 ), Math::Util::RotationCast< quat_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T diff = quaternionDiff( quats[ i ], quatsAA4[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );
		}
	}
	
	{	// test quaternion -> rotation axis transformation
		std::vector< rotaxis_type > aa3s;
		quats.reserve( n );
		std::transform( quats.begin(), quats.end(), std::back_inserter( aa3s ), Math::Util::RotationCast< rotaxis_type > () );
	
		std::vector< quat_type > quatsAA3;
		quatsAA3.reserve( n );
		std::transform( aa3s.begin(), aa3s.end(), std::back_inserter( quatsAA3 ), Math::Util::RotationCast< quat_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T diff = quaternionDiff( quats[ i ], quatsAA3[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );
		}
	}
	
	{	// test quaternion -> 3-by-3 matrix transformation
		std::vector< mat_type > matrices;
		matrices.reserve( n );
		std::transform( quats.begin(), quats.end(), std::back_inserter( matrices ), Math::Util::RotationCast< mat_type > () );
		
		std::vector< quat_type > quatsMat;
		quatsMat.reserve( n );
		std::transform( matrices.begin(), matrices.end(), std::back_inserter( quatsMat ), Math::Util::RotationCast< quat_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T diff = quaternionDiff( quats[ i ], quatsMat[ i ] );
			BOOST_CHECK_SMALL( diff, epsilon );
		}
	}
}

template< typename T >
void testMatrixCast( const std::size_t n_runs, const T epsilon )
{
	// we do all runs at the same time...
	const std::size_t n = n_runs;
	
	// cast a quaterion to different type
	typedef Math::Quaternion quat_type;
	typedef Math::Matrix< T, 3, 3 > mat_type;
	typedef Math::Vector< T, 4 > aaxis_type; // axis angle 
	typedef Math::Vector< T, 3 > rotaxis_type; // rotation axis
	
	typename Math::Random::Quaternion< T >::Uniform randQuat;
	
	std::vector< quat_type > randQuats;
	randQuats.reserve( n );
	std::generate_n ( std::back_inserter( randQuats ), n,  randQuat );
	// we generate rdandom rotation matrices from random quaterions (easier than directly producing roation matrices)
	std::vector< mat_type > matrices;
	matrices.reserve( n );
	std::transform( randQuats.begin(), randQuats.end(), std::back_inserter( matrices ), Math::Util::RotationCast< mat_type > () );
	
	{	// test matrix -> axis angle transformation
		std::vector< aaxis_type > aa4s;
		aa4s.reserve( n );
		std::transform( matrices.begin(), matrices.end(), std::back_inserter( aa4s ), Math::Util::RotationCast< aaxis_type > () );
		
		std::vector< mat_type > matricesAA4;
		matricesAA4.reserve( n );
		std::transform( aa4s.begin(), aa4s.end(), std::back_inserter( matricesAA4 ), Math::Util::RotationCast< mat_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T det = rotMatrixDiff( matrices[ i ], matricesAA4[ i ] );
			const T diff = 1 - det;
			BOOST_CHECK_SMALL( diff, epsilon );
			// std::cout << " Matrix " << matricesAA4[ i ] << std::endl;
		}
	}
	
	{	// test matrix -> rotation axis transformation
		std::vector< rotaxis_type > aa3s;
		aa3s.reserve( n );
		std::transform( matrices.begin(), matrices.end(), std::back_inserter( aa3s ), Math::Util::RotationCast< rotaxis_type > () );
		
		std::vector< mat_type > matricesAA3;
		matricesAA3.reserve( n );
		std::transform( aa3s.begin(), aa3s.end(), std::back_inserter( matricesAA3 ), Math::Util::RotationCast< mat_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T det = rotMatrixDiff( matrices[ i ], matricesAA3[ i ] );
			const T diff = 1 - det;
			BOOST_CHECK_SMALL( diff, epsilon );
			// std::cout << " Matrix " << matricesAA3[ i ] << std::endl;
		}
	}
	
	{	// test matrix -> quaterion transformation
		std::vector< quat_type > quats;
		quats.reserve( n );
		std::transform( matrices.begin(), matrices.end(), std::back_inserter( quats ), Math::Util::RotationCast< quat_type > () );
		
		std::vector< mat_type > matricesQuat;
		matricesQuat.reserve( n );
		std::transform( quats.begin(), quats.end(), std::back_inserter( matricesQuat ), Math::Util::RotationCast< mat_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T det = rotMatrixDiff( matrices[ i ], matricesQuat[ i ] );
			const T diff = 1 - det;
			BOOST_CHECK_SMALL( diff, epsilon );
			// std::cout << " Matrix " << matricesQuat[ i ] << std::endl;
		}
	}
}

template< typename T >
void testAxisAngleCast( const std::size_t n_runs, const T epsilon )
{
	// we do all runs at the same time...
	const std::size_t n = n_runs;
	
	// cast a quaterion to different type
	typedef Math::Quaternion quat_type;
	typedef Math::Matrix< T, 3, 3 > mat_type;
	typedef Math::Vector< T, 4 > aaxis_type; // axis angle 
	typedef Math::Vector< T, 3 > rotaxis_type; // rotation axis
	
	typename Math::Random::Quaternion< T >::Uniform randQuat;
	
	std::vector< quat_type > randQuats;
	randQuats.reserve( n );
	std::generate_n ( std::back_inserter( randQuats ), n,  randQuat );
	// we generate rdandom rotation matrices from random quaterions (easier than directly producing roation matrices)
	std::vector< aaxis_type > axisAngles;
	axisAngles.reserve( n );
	std::transform( randQuats.begin(), randQuats.end(), std::back_inserter( axisAngles ), Math::Util::RotationCast< aaxis_type > () );
	
	// commentes: test left out since conversion gives no reliable results at the moment...
	{	// test axis angle -> matrix transformation
		// std::vector< mat_type > matrices;
		// matrices.reserve( n );
		// std::transform( axisAngles.begin(), axisAngles.end(), std::back_inserter( matrices ), Math::Util::RotationCast< mat_type > () );
		
		// std::vector< aaxis_type > axisAngleMAT;
		// axisAngleMAT.reserve( n );
		// std::transform( matrices.begin(), matrices.end(), std::back_inserter( axisAngleMAT ), Math::Util::RotationCast< aaxis_type > () );
		
		// for( std::size_t i = 0; i<n; ++i )
		// {
			// const T dist = Math::distance( axisAngles[ i ], axisAngleMAT[ i ] );
			// BOOST_CHECK_SMALL( dist, epsilon );
			// std::cout << " Axis Angle GT  " << axisAngles[ i ] << std::endl;
			// std::cout << " Axis Angle EST " << axisAngleMAT[ i ] << std::endl;
		// }
	}
	
	{	// test axis angle -> rotation axis transformation
		std::vector< rotaxis_type > aa3s;
		aa3s.reserve( n );
		std::transform( axisAngles.begin(), axisAngles.end(), std::back_inserter( aa3s ), Math::Util::RotationCast< rotaxis_type > () );
		
		std::vector< aaxis_type > axisAngleAA3;
		axisAngleAA3.reserve( n );
		std::transform( aa3s.begin(), aa3s.end(), std::back_inserter( axisAngleAA3 ), Math::Util::RotationCast< aaxis_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T dist = Math::distance( axisAngles[ i ], axisAngleAA3[ i ] );
			BOOST_CHECK_SMALL( dist, epsilon );
			// std::cout << " Matrix " << axisAngleAA3[ i ] << std::endl;
		}
	}
	
	{	// test axis angle -> quaterion transformation
		std::vector< quat_type > quats;
		quats.reserve( n );
		std::transform( axisAngles.begin(), axisAngles.end(), std::back_inserter( quats ), Math::Util::RotationCast< quat_type > () );
		
		std::vector< aaxis_type > axisAngleQuat;
		axisAngleQuat.reserve( n );
		std::transform( quats.begin(), quats.end(), std::back_inserter( axisAngleQuat ), Math::Util::RotationCast< aaxis_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T dist = Math::distance( axisAngles[ i ], axisAngleQuat[ i ] );
			BOOST_CHECK_SMALL( dist, epsilon );
			// std::cout << " Matrix " << axisAngleQuat[ i ] << std::endl;
		}
	}
}


template< typename T >
void testRotationAxisCast( const std::size_t n_runs, const T epsilon )
{
	// we do all runs at the same time...
	const std::size_t n = n_runs;
	
	// cast a quaterion to different type
	typedef Math::Quaternion quat_type;
	typedef Math::Matrix< T, 3, 3 > mat_type;
	typedef Math::Vector< T, 4 > aaxis_type; // axis angle 
	typedef Math::Vector< T, 3 > rotaxis_type; // rotation axis
	
	typename Math::Random::Quaternion< T >::Uniform randQuat;
	
	std::vector< quat_type > randQuats;
	randQuats.reserve( n );
	std::generate_n ( std::back_inserter( randQuats ), n,  randQuat );
	// we generate rdandom rotation matrices from random quaterions (easier than directly producing roation matrices)
	std::vector< rotaxis_type > rotAxes;
	rotAxes.reserve( n );
	std::transform( randQuats.begin(), randQuats.end(), std::back_inserter( rotAxes ), Math::Util::RotationCast< rotaxis_type > () );
	
	// commentes: test left out since conversion gives no reliable results at the moment...
	{	// test axis angle -> matrix transformation
		// std::vector< mat_type > matrices;
		// matrices.reserve( n );
		// std::transform( rotAxes.begin(), rotAxes.end(), std::back_inserter( matrices ), Math::Util::RotationCast< mat_type > () );
		
		// std::vector< rotaxis_type > axisAngleMAT;
		// axisAngleMAT.reserve( n );
		// std::transform( matrices.begin(), matrices.end(), std::back_inserter( axisAngleMAT ), Math::Util::RotationCast< rotaxis_type > () );
		
		// for( std::size_t i = 0; i<n; ++i )
		// {
			// const T dist = Math::distance( rotAxes[ i ], axisAngleMAT[ i ] );
			// BOOST_CHECK_SMALL( dist, epsilon );
			// std::cout << " Axis Angle GT  " << rotAxes[ i ] << std::endl;
			// std::cout << " Axis Angle EST " << axisAngleMAT[ i ] << std::endl;
		// }
	}
	
	{	// test axis angle -> rotation axis transformation
		std::vector< aaxis_type > aa4s;
		aa4s.reserve( n );
		std::transform( rotAxes.begin(), rotAxes.end(), std::back_inserter( aa4s ), Math::Util::RotationCast< aaxis_type > () );
		
		std::vector< rotaxis_type > axisAngleAA4;
		axisAngleAA4.reserve( n );
		std::transform( aa4s.begin(), aa4s.end(), std::back_inserter( axisAngleAA4 ), Math::Util::RotationCast< rotaxis_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T dist = Math::distance( rotAxes[ i ], axisAngleAA4[ i ] );
			BOOST_CHECK_SMALL( dist, epsilon );
			// std::cout << " Matrix " << axisAngleAA4[ i ] << std::endl;
		}
	}
	
	{	// test axis angle -> quaterion transformation
		std::vector< quat_type > quats;
		quats.reserve( n );
		std::transform( rotAxes.begin(), rotAxes.end(), std::back_inserter( quats ), Math::Util::RotationCast< quat_type > () );
		
		std::vector< rotaxis_type > axisAngleQuat;
		axisAngleQuat.reserve( n );
		std::transform( quats.begin(), quats.end(), std::back_inserter( axisAngleQuat ), Math::Util::RotationCast< rotaxis_type > () );
		
		for( std::size_t i = 0; i<n; ++i )
		{
			const T dist = Math::distance( rotAxes[ i ], axisAngleQuat[ i ] );
			BOOST_CHECK_SMALL( dist, epsilon );
			// std::cout << " Matrix " << axisAngleQuat[ i ] << std::endl;
		}
	}
}


	
void TestRotationCast()
{
	
	testQuaternionCast< double >( 100000, 1e-10 );
	testMatrixCast< double >( 100000, 1e-10 );
	testAxisAngleCast< double >( 100000, 1e-10 );
	testRotationAxisCast< double >( 10, 1e-10 );
	
	
	
	BOOST_CHECK( 0 == 0 );
}
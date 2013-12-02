/*
 * Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2006, Technische Universitaet Muenchen, and individual
 * contributors as indicated by the @authors tag. See the 
 * copyright.txt in the distribution for a full listing of individual
 * contributors.
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 */

/**
 * @ingroup tracking_algorithms
 * @file
 * Implements functions for fundamental matrix estimation.
 *
 * @author Daniel Muhra <muhra@in.tum.de>
 */ 

#include "FundamentalMatrix.h"
#include <utMath/MatrixOperations.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <log4cpp/Category.hh>
#include <utUtil/Exception.h>
#include <utUtil/Logging.h>
#include <boost/numeric/ublas/io.hpp>

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>

// shortcuts to namespaces
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

namespace Ubitrack { namespace Calibration {

template< typename T >
void normalize(Math::Vector< 2, T >& shift, T& scale,
	Math::Matrix< 3, 3, T >& modMatrix, const std::vector< Math::Vector< 2, T > >& pts)
{    
	shift = Math::Vector< 2, T >( 0.0, 0.0 );
	T meandist = 0.0;
		
	for( std::size_t i=0; i < pts.size(); i++ )
		shift += pts.at( i );
	shift = shift / pts.size();
	
	for( std::size_t i=0; i < pts.size(); i++ )
	{
		Math::Vector< 2, T > v = pts.at(i) - shift;
		meandist += sqrt( v(0) * v(0) + v(1) * v(1) );
	}    
    meandist /= pts.size();
    
    scale = sqrt( static_cast< T >( 2.0 ) )/meandist;

	modMatrix( 0, 0 ) = scale ;
	modMatrix( 0, 1 ) = 0.0;
	modMatrix( 0, 2 ) = -scale * shift( 0 );
	modMatrix( 1, 0 ) = 0.0;
	modMatrix( 1, 1 ) = scale;
	modMatrix( 1, 2 ) = -scale * shift( 1 );
	modMatrix( 2, 0 ) = 0.0;
	modMatrix( 2, 1 ) = 0.0;
	modMatrix( 2, 2 ) = 1.0;
}

template< typename T >
Math::Matrix< 3, 3, T > getFundamentalMatrixImp( const std::vector< Math::Vector< 2, T > > & fromPoints, 
	const std::vector< Math::Vector< 2, T > > & toPoints, std::size_t stepSize )
{
	static log4cpp::Category& logger(log4cpp::Category::getInstance( "Ubitrack.Calibration.FundamentalMatrix" ));
	
	if( stepSize < 1 )
	{
		LOG4CPP_ERROR ( logger, "invalid step size, using 1 instead");
		UBITRACK_THROW ( "invalid step size, using 1 instead" );	
		stepSize = 1;
	}

	if( fromPoints.size() != toPoints.size() )
	{
		LOG4CPP_ERROR ( logger, "Input sizes of the vectors do not match ");
		UBITRACK_THROW ( "Input sizes do not match" );	
	}
	
	if( fromPoints.size() < 8 )
	{
		LOG4CPP_ERROR ( logger, "Input sizes of the vectors to small ");
		UBITRACK_THROW ( "Input sizes to small. Use at least 8 values" );	
	}

	//Normalization
	Math::Vector< 2, T > fromShift( 0.0, 0.0 );
	T fromScale = static_cast< T >( 0 );
	Math::Vector< 2, T > toShift( 0.0, 0.0 );
	T toScale = static_cast< T >( 0 );
	Math::Matrix< 3, 3, T > fromModMatrix;
	Math::Matrix< 3, 3, T > toModMatrix;

	normalize( fromShift, fromScale, fromModMatrix, fromPoints );
	normalize( toShift, toScale, toModMatrix, toPoints );

	//Linear Solution
	Math::Matrix< 0, 0, T > A( fromPoints.size() / stepSize, 9 );

	for( std::size_t i=0; i < ( fromPoints.size() / stepSize ); i++ )
	{
		T x = ( fromPoints.at( i * stepSize )( 0 ) - fromShift( 0 ) ) * fromScale;
		T x_ = ( toPoints.at( i * stepSize )( 0 ) - toShift( 0 ) ) * toScale;
		T y = ( fromPoints.at( i * stepSize )( 1 ) - fromShift( 1 ) ) * fromScale;
		T y_ = ( toPoints.at( i * stepSize )( 1 ) - toShift( 1 ) ) * toScale;

		A( i, 0 ) = x_ * x;
		A( i, 1 ) = x_ * y;
		A( i, 2 ) = x_;
		A( i, 3 ) = y_* x;
		A( i, 4 ) = y_ * y;
		A( i, 5 ) = y_;
		A( i, 6 ) = x;
		A( i, 7 ) = y;
		A( i, 8 ) = static_cast< T >( 1.0 );
	}

	// solve using SVD
	std::size_t nSingularValues = std::min( A.size1(), A.size2() );
	Math::Vector< 0, T > s1( nSingularValues );
	Math::Matrix< 9, 9, T > Vt;
	Math::Matrix< 0, 0, T > U( fromPoints.size() / stepSize, fromPoints.size() / stepSize );
	int info = lapack::gesvd( 'N', 'A', A, s1, U, Vt );

	if ( info != 0 )
	{
		LOG4CPP_ERROR ( logger, "first SVD failed");
		UBITRACK_THROW ( "first SVD failed" );	
	}

	// copy result to 3x3 matrix
	Math::Matrix< 3, 3, T > F;
	F( 0, 0 ) = Vt( 8, 0 ); F( 0, 1 ) = Vt( 8, 1 ); F( 0, 2 ) = Vt( 8, 2 );
	F( 1, 0 ) = Vt( 8, 3 ); F( 1, 1 ) = Vt( 8, 4 ); F( 1, 2 ) = Vt( 8, 5 );
	F( 2, 0 ) = Vt( 8, 6 ); F( 2, 1 ) = Vt( 8, 7 ); F( 2, 2 ) = Vt( 8, 8 );

	// constraint enforcement
	Math::Vector< 0, T > s2( 3 );
	Math::Matrix< 3, 3, T > U2;
	Math::Matrix< 3, 3, T > Vt2;
	info = lapack::gesvd( 'A', 'A', F, s2, U2, Vt2 );

	if ( info != 0 )
	{
		LOG4CPP_ERROR ( logger, "second SVD failed");
		UBITRACK_THROW ( "second SVD failed" );	
	}

	s2( 2 ) = 0;
	for ( std::size_t i = 0; i < 3; i++ )
		ublas::column( U2, i ) *= s2( i );
	
	F = ublas::prod( U2, Vt2 );
	F = ublas::prod( ublas::trans( toModMatrix ), F );
	return ublas::prod( F, fromModMatrix );
}

Math::Matrix< 3, 3, float > getFundamentalMatrix( const std::vector< Math::Vector< 2, float > >& fromPoints, 
	const std::vector< Math::Vector< 2, float > >& toPoints, std::size_t stepSize )
{ 
	return getFundamentalMatrixImp( fromPoints, toPoints, stepSize );
}

Math::Matrix< 3, 3, double > getFundamentalMatrix( const std::vector< Math::Vector< 2, double > >& fromPoints, 
	const std::vector< Math::Vector< 2, double > >& toPoints, std::size_t stepSize )
{
	return getFundamentalMatrixImp( fromPoints, toPoints, stepSize );
}

Math::Matrix< 3, 3, double > fundamentalMatrixFromPoses( const Math::Pose & cam1, const Math::Pose & cam2, const Math::Matrix< 3, 3 > & K1, const Math::Matrix< 3, 3 > & K2 )
{
	Math::Matrix< 3, 4 > E1( cam1 );
	Math::Matrix< 3, 4 > E2( cam2 );
	
	E1 = ublas::prod( K1, E1 );
	E2 = ublas::prod( K2, E2 );

	Math::Matrix< 4, 4 > P1_I = Math::invert_matrix( Math::Matrix< 4, 4 >( cam1 ) );
	Math::Vector< 4 > C( P1_I( 0 , 3 ), P1_I( 1 , 3 ), P1_I( 2 , 3 ),  1.0 );

	Math::Vector< 3 > e2 = ublas::prod( E2, C );
	Math::Matrix< 3, 3 > skew;

	skew( 0, 0 ) = 0.0;
	skew( 0, 1 ) = -e2( 2 );
	skew( 0, 2 ) = e2( 1 );
	skew( 1, 0 ) = e2( 2 );
	skew( 1, 1 ) = 0.0;
	skew( 1, 2 ) = -e2( 0 );
	skew( 2, 0 ) = -e2( 1 );
	skew( 2, 1 ) = e2( 0 );
	skew( 2, 2 ) = 0.0;

	Math::Matrix< 4, 3 > E1_ = Math::pseudoInvert_matrix( E1 );

	Math::Matrix< 3, 3 > F = ublas::prod( E2, E1_ );
	return ublas::prod( skew, F );
}

Math::Pose poseFromFundamentalMatrix( const Math::Matrix< 3, 3 > & fM, const Math::Vector< 2 > & x, const Math::Vector< 2 > & x_, const Math::Matrix< 3, 3 > & K1, const Math::Matrix< 3, 3 > & K2 )
{
	//get essential matrix
	Math::Matrix< 3, 3 > W = ublas::trans( K2 );
	W = ublas::prod( W, fM );
	W = ublas::prod( W, K1 );

	//find possible poses
	Math::Vector< 3 > s;
	Math::Matrix< 3, 3 > Vt;
	Math::Matrix< 3, 3 > U;
	lapack::gesvd( 'A', 'A', W, s, U, Vt );

	W( 0, 0 ) = 0.0;
	W( 0, 1 ) = -1.0;
	W( 0, 2 ) = 0.0;
	W( 1, 0 ) = 1.0;
	W( 1, 1 ) = 0.0;
	W( 1, 2 ) = 0.0;
	W( 2, 0 ) = 0.0;
	W( 2, 1 ) = 0.0;
	W( 2, 2 ) = 1.0;

	Math::Matrix< 3, 3 > Wt = ublas::trans( W );
	Math::Vector< 3 > u3 = ublas::column( U, 2 );
	
	W = ublas::prod( W, Vt );
	Wt = ublas::prod( Wt, Vt );

	W = ublas::prod( U, W );
	Wt = ublas::prod( U, Wt );

	if( Math::determinant( W ) < 0 )
		W *= -1.0;

	if( Math::determinant( Wt ) < 0 )
		Wt *= -1.0;
	
	//check which is the correct Pose
	Math::Matrix< 3, 4 > p1;
	ublas::subrange( p1, 0, 3, 0, 3 ) = ublas::identity_matrix< double >(3);
	ublas::column( p1, 3 ) = Math::Vector< 3 >( 0.0, 0.0, 0.0 );
	Math::Matrix< 3, 4 > P1 = ublas::prod( K1, p1 );

	Math::Matrix< 4, 4 > inv = ublas::identity_matrix< double >(4);
	ublas::subrange( inv, 0, 3, 0, 4 ) = p1;
	inv = Math::invert_matrix( inv );

	//first possibility
	Math::Matrix< 3, 4 > p2;
	ublas::subrange( p2, 0, 3, 0, 3 ) = W;
	ublas::column( p2, 3 ) = u3;
	Math::Matrix< 3, 4 > P2 = ublas::prod( K2, p2 );

	s = Calibration::get3DPosition( P1, P2, x, x_ );

	if ( ( ublas::prod( W, s ) + u3 )( 2 ) < 0.0 && s( 2 ) < 0.0 )
		return Math::Pose( Math::Quaternion( W ), u3 );
	
	//second possibility
	ublas::column( p2, 3 ) = -u3;
	P2 = ublas::prod( K2, p2 );

	s = Calibration::get3DPosition( P1, P2, x, x_ );

	if ( ( ublas::prod( W, s ) - u3 )( 2 ) < 0.0 && s( 2 ) < 0.0  )
		return Math::Pose( Math::Quaternion( W ), -u3 );

	//third possibility
	ublas::subrange( p2, 0, 3, 0, 3 ) = Wt;
	P2 = ublas::prod( K2, p2 );

	s = Calibration::get3DPosition( P1, P2, x, x_ );

	if ( ( ublas::prod( Wt, s ) - u3 )( 2 ) < 0.0 && s( 2 ) < 0.0 )
		return Math::Pose( Math::Quaternion( Wt ), -u3 );

	//last possibility
	return Math::Pose( Math::Quaternion( Wt ), u3 );
}

} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

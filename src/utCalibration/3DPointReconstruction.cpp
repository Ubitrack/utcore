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
 * Implements  functions for 3D Point Reconstruction
 * @author Daniel Muhra <muhra@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#include "3DPointReconstruction.h"
#include <iostream>
#include <iterator>

#include <utUtil/Logging.h>
#include <utUtil/Exception.h>
#include <utMath/GaussNewton.h>
#include <utMath/LevenbergMarquardt.h>
#include <utCalibration/Function/SinglePointMultiProjection.h>


#include <log4cpp/Category.hh>

namespace ublas = boost::numeric::ublas;

#ifdef HAVE_LAPACK
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <utCalibration/Munkres.h>
namespace lapack = boost::numeric::bindings::lapack;
#endif

namespace Ubitrack { namespace Calibration {

/** internal of pointToPointDist function */
template< typename T >
T pointToPointDistImp( const Math::Vector< 3, T > & from, const Math::Vector< 3, T > & to, const Math::Matrix< 3, 3, T > & fM  )
{
	Math::Vector< 3, T > from_  = ublas::prod( fM, from );

	T term = from_( 0 ) * to( 0 ) + from_( 1 )* to( 1 ) + from_( 2 );
	return ( term * term ) / ( from_( 0 ) * from_( 0 ) + from_( 1 )* from_( 1 ) );
}

float pointToPointDist( const Math::Vector< 2, float > & from, const Math::Vector< 2, float > & to, const Math::Matrix< 3, 3, float > & fM  )
{
	const Math::Vector< 3, float > from_( from( 0 ), from( 1 ), 1.0f );
	const Math::Vector< 3, float > to_( to( 0 ), to( 1 ), 1.0f );

	return pointToPointDistImp( from_, to_, fM );
}

double pointToPointDist( const Math::Vector< 2, double > & from, const Math::Vector< 2, double > & to, const Math::Matrix< 3, 3, double > & fM  )
{
	const Math::Vector< 3, double > from_( from( 0 ), from( 1 ), 1.0 );
	const Math::Vector< 3, double > to_( to( 0 ), to( 1 ), 1.0 );

	return pointToPointDistImp( from_, to_, fM );
}

float pointToPointDist( const Math::Vector< 3, float > & from, const Math::Vector< 3, float > & to, const Math::Matrix< 3, 3, float > & fM  )
{
	return pointToPointDistImp( from, to, fM );
}

double pointToPointDist( const Math::Vector< 3, double > & from, const Math::Vector< 3, double > & to, const Math::Matrix< 3, 3, double > & fM  )
{
	return pointToPointDistImp( from, to, fM );
}

/** internal of get3DPostion function */
#ifdef HAVE_LAPACK
template< typename T >
Math::Vector< 3, T > get3DPositionImp( const Math::Matrix< 3, 4, T > & P1, const Math::Matrix< 3, 4, T > & P2, const Math::Vector< 2, T > & x, const Math::Vector< 2, T > & x_ )
{
	Math::Matrix< 4, 4, T > A;
	ublas::row( A, 0 ) = ( x( 0 ) * ublas::row( P1, 2 ) ) - ublas::row( P1, 0 );
	ublas::row( A, 1 ) = ( x( 1 ) * ublas::row( P1, 2 ) ) - ublas::row( P1, 1 );
	ublas::row( A, 2 ) = ( x_( 0 ) * ublas::row( P2, 2 ) ) - ublas::row( P2, 0 );
	ublas::row( A, 3 ) = ( x_( 1 ) * ublas::row( P2, 2 ) ) - ublas::row( P2, 1 );

	//solving using svd
	ublas::vector< T > s1( 4 );
	Math::Matrix< 4, 4, T > Vt;
	Math::Matrix< 4, 4, T > U;
	lapack::gesvd( 'N', 'A', A, s1, U, Vt );

	return ( Math::Vector< 3, T >( Vt( 3, 0 ) / Vt( 3, 3), Vt( 3, 1 ) / Vt( 3, 3), Vt( 3, 2 ) / Vt( 3, 3) ) );
}

Math::Vector< 3, float > get3DPosition( const Math::Matrix< 3, 4, float > & P1, const Math::Matrix< 3, 4, float > & P2, const Math::Vector< 2, float > & x, const Math::Vector< 2, float > & x_ )
{
	return get3DPositionImp( P1, P2, x, x_ );
}

Math::Vector< 3, double > get3DPosition( const Math::Matrix< 3, 4, double > & P1, const Math::Matrix< 3, 4, double > & P2, const Math::Vector< 2, double > & x, const Math::Vector< 2, double > & x_ )
{
	return get3DPositionImp( P1, P2, x, x_ );
}

/** internal of reconstruct3DPoints function */
template< typename ForwardIterator1, typename ForwardIterator2 >
Math::Vector< 3, typename std::iterator_traits< ForwardIterator2 >::value_type::value_type  > get3DPositionImpl( const ForwardIterator1 iBegin, const ForwardIterator1 iEnd, ForwardIterator2 iPoints )
{
	// setting shortcut to Type= double or float, kind of ugly to read
	// since it is a value_type of a value_type :-(
	typedef typename std::iterator_traits< ForwardIterator2 >::value_type::value_type Type;

	//determening the size
	unsigned n ( iEnd - iBegin );
	if( n < 2 )
		UBITRACK_THROW( "3d point estimation requires at least 2 matrices and 2 image points." );

	ublas::matrix< Type, ublas::column_major > A( n * 3, 4 );
	
	unsigned i( 0 );
	for ( ForwardIterator1 it ( iBegin ); it != iEnd; ++i, ++it, ++iPoints )
	{
		// building Matrix to solve for null space
		// each 3x4 block i contains information of the form
		// A_i = [ x_i, y_i, 1 ]_x * P_i
		Math::Matrix< 3, 3, Type > skew;
		ublas::row( skew, 0 ) = Math::Vector< 3, Type >( 0, 1, -((*iPoints)( 1 )) );
		ublas::row( skew, 1 ) = Math::Vector< 3, Type >( -1, 0, (*iPoints)( 0 ) );
		ublas::row( skew, 2 ) = Math::Vector< 3, Type >( (*iPoints)( 1 ), -((*iPoints)( 0 )), 0 );
		ublas::subrange( A, i*3, i*3+3, 0, 4 ) = ublas::prod( skew, (*it) );
	}

	ublas::vector< Type > s( 4 );
	Math::Matrix< 4, 4, Type > Vt;
	ublas::matrix< Type, ublas::column_major > U( 3*n, 3*n );
	if( lapack::gesvd( 'N', 'A', A, s, U, Vt ) != 0 )
		UBITRACK_THROW ( "SVD for point reconstruction failed." );
		
	Math::Vector< 4, Type > vec = ublas::row( Vt, 3 );
	
	for ( ForwardIterator1 it ( iBegin ); it != iEnd; it++ )
	{
		// check for the direction of the point
		Type s = ( (*it)( 2, 0 )*vec( 0 ) + (*it)( 2, 1 )*vec( 1 ) + (*it)( 2, 2 )*vec( 2 ) + (*it)( 2, 3 )*vec( 3 ) );
		 if( 0 > s )
		{
			vec = -vec;
			break;
		}	
	}
	vec /= vec( 3 );
	return Math::Vector< 3, Type >( ublas::subrange( vec, 0, 3 ) );
}

/** internal non-linear optimization function for 3d point reconstruction */
template< typename ForwardIterator1, typename ForwardIterator2 >
Math::Vector< 3, typename std::iterator_traits< ForwardIterator1 >::value_type::value_type  > optimize3DPositionImpl( const ForwardIterator1 iBegin, const ForwardIterator1 iEnd, ForwardIterator2 iPoints, const Math::Vector< 3, typename std::iterator_traits< ForwardIterator2 >::value_type::value_type > &initialPoint )
{
	// shortcut to double/float
	typedef typename std::iterator_traits< ForwardIterator1 >::value_type::value_type Type;
	
	unsigned n ( iEnd - iBegin );
	Function::SinglePointMultiProjection< Type, ForwardIterator1 > func( iBegin, iEnd );
	
	// prepare the image measurement vector for the minimization
	ublas::vector< Type > measurement( n * 2 );
	ForwardIterator2 it( iPoints );
	for ( unsigned i ( 0 ); i < n; ++i, ++it )
		ublas::subrange( measurement, i*2, (i*2)+2 ) = *it;

	// prepare the input 3-vector to be optimized
	ublas::vector< Type > parameters( 3 );
	parameters = initialPoint;
	
	// perform optimization
	Type residual = Ubitrack::Math::levenbergMarquardt( func, parameters, measurement, Math::OptTerminate( 200, 1e-6 ), Math::OptNoNormalize() );

	return Math::Vector< 3 , Type >( parameters );
	
}

Math::Vector< 3, float > get3DPosition( const std::vector< Math::Matrix< 3, 4, float > > &P, const std::vector< Math::Vector< 2, float > > &points, unsigned flag )
{
	if( P.size() != points.size() )
		UBITRACK_THROW( "no equal amount of camera projections and corresponding points." );
	Math::Vector< 3, float > result = get3DPositionImpl( P.begin(), P.end(), points.begin() );
	if( flag > 0 )
		result = optimize3DPositionImpl( P.begin(), P.end(), points.begin(), result );
	return result;
}

Math::Vector< 3, double > get3DPosition( const std::vector< Math::Matrix< 3, 4, double > > &P, const std::vector< Math::Vector< 2, double > >& points, unsigned flag )
{
	if( P.size() != points.size() )
		UBITRACK_THROW( "no equal amount of camera projections and corresponding points." );
	Math::Vector< 3, double > result = get3DPositionImpl( P.begin(), P.end(), points.begin() );
	if( flag > 0 )
		result = optimize3DPositionImpl( P.begin(), P.end(), points.begin(), result );
	
	return result;
}


/** internal of reconstruct3DPoints function */
template< typename T >
std::vector< Math::Vector< 3, T > > reconstruct3DPointsImp( const std::vector< Math::Vector< 2, T > > & p1, const std::vector< Math::Vector< 2, T > > & p2,
																			const Math::Matrix< 3, 4, T > & P1, const Math::Matrix< 3, 4, T > & P2, const Math::Matrix< 3, 3, T > & fM )
{
	
	unsigned int p1Size = p1.size();
	unsigned int p2Size = p2.size();

	//create match matrix
	ublas::matrix< T > matrix( p1Size, p2Size );

	for( unsigned int row=0; row < p1Size; row++ )
	{
		for( unsigned int col=0; col < p2Size; col++ )
		{
			matrix( row, col ) = pointToPointDist( p1.at( row ), p2.at( col ), fM );
		}
	}

	Munkres< T > m( matrix );
	m.solve();
	std::vector< unsigned int > matchList = m.getRowMatchList();

	std::vector< Math::Vector< 3, T > > list;

	for(unsigned int i=0; i < p1Size; i++ )
	{
		if( matchList.at( i ) < p2Size )
		{
			list.push_back( get3DPosition( P1, P2, p1.at( i ), p2.at( matchList.at( i ) ) ) );
		}
	}

	return list;
}

std::vector< Math::Vector< 3, float > > reconstruct3DPoints( const std::vector< Math::Vector< 2, float > > & p1, const std::vector< Math::Vector< 2, float > > & p2,
																			const Math::Matrix< 3, 4, float > & P1, const Math::Matrix< 3, 4, float > & P2, const Math::Matrix< 3, 3, float > & fM )
{
	return reconstruct3DPointsImp( p1, p2, P1, P2, fM );
}

std::vector< Math::Vector< 3, double > > reconstruct3DPoints( const std::vector< Math::Vector< 2, double > > & p1, const std::vector< Math::Vector< 2, double > > & p2,
																			const Math::Matrix< 3, 4, double > & P1, const Math::Matrix< 3, 4, double > & P2, const Math::Matrix< 3, 3, double > & fM )
{
	return reconstruct3DPointsImp( p1, p2, P1, P2, fM );
}
#endif

} } // namespace Ubitrack::Calibration

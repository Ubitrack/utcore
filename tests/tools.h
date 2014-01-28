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


#ifndef __TESTS_TOOLS_H_INCLUDED__
#define __TESTS_TOOLS_H_INCLUDED__

#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utMath/Quaternion.h>
#include <utMath/Blas1.h>

#include <math.h>
#include <numeric> // std::accumulate

template< class T > 
T random( T a, T b )
{ return T( rand() ) / T( RAND_MAX ) * ( b - a ) + a; }

template< class M > 
void randomMatrix( M& m )
{
	for ( unsigned i = 0; i < m.size1(); i++ )
		for ( unsigned j = 0; j < m.size2(); j++ )
			m( i, j ) = random( typename M::value_type( -100 ), typename M::value_type( 100 ) );
}

template< typename T, std::size_t N > 
Ubitrack::Math::Vector< T, N > randomVector( T maxVal = 100.0 )
{
	Ubitrack::Math::Vector< T, N > v;
	for ( std::size_t i = 0; i < v.size(); i++ )
		v( i ) = random( T( -maxVal ), T( maxVal ) );
	return v;
}

template< class MA, class MB > 
typename MA::value_type matrixDiff( const MA& ma, const MB& mb )
{
	typedef typename MA::value_type return_type;
	return_type d = 0.0;
	for ( std::size_t i = 0; i < ma.size1(); i++ )
		for ( std::size_t j = 0; j < ma.size2(); j++ )
			d += fabs( ma( i, j ) - mb( i, j ) );
	return d / boost::numeric::ublas::norm_frobenius( ma );
}

template< class VA, class VB > 
typename VA::value_type vectorDiff( const VA& va, const VB& vb )
{
	typedef typename VA::value_type return_type;
	return_type d = 0.0;
	for ( std::size_t i = 0; i < va.size(); i++ )
		d += fabs( va( i ) - vb( i ) );
	return d / boost::numeric::ublas::norm_2( va );
}

template< class VA, class VB > 
typename VA::value_type vectorDistance( const VA& va, const VB& vb )
{
	typedef typename VA::value_type return_type;
	return_type d = 0.0;
	for ( std::size_t i = 0; i < va.size(); i++ )
		d += fabs( va( i ) - vb( i ) );
	return std::sqrt( d );
}

template< class MA, class MB > 
typename MA::value_type homMatrixDiff( const MA& A, const MB& B )
{

	typedef typename MA::value_type value_type;
	// normalize both H and Htest
	// ublas is evil. Don't try this: Htest /= ublas::norm_frobenius( Htest )
	const value_type normA = boost::numeric::ublas::norm_frobenius( A );
	const value_type normB = boost::numeric::ublas::norm_frobenius( B );
	
	value_type dp = 0.0; // sum of differences A-B
	value_type dm = 0.0; // sum of differences A+B (in case A ~ -B)
	for ( std::size_t i = 0; i < A.size1(); i++ )
		for ( std::size_t j = 0; j < A.size2(); j++ )
		{
			dp += std::fabs( A( i, j ) / normA - B( i, j ) / normB );
			dm += std::fabs( A( i, j ) / normA + B( i, j ) / normB );
		}
		
	return std::min( dp, dm );
}

static Ubitrack::Math::Quaternion randomQuaternion()
{ 
	return Ubitrack::Math::Quaternion( random( -1.0, 1.0 ), random( -1.0, 1.0 ), 
		random( -1.0, 1.0 ), random( -1.0, 1.0 ) ).normalize(); 
}

static double quaternionDiff( const Ubitrack::Math::Quaternion& _a, const Ubitrack::Math::Quaternion& _b )
{
	Ubitrack::Math::Quaternion a( _a ); a.normalize();
	Ubitrack::Math::Quaternion b( _b ); b.normalize();
	if ( a.x() * b.x() + a.y() * b.y() + a.z() * b.z() + a.w() * b.w() < 0 )
		a = -a;
	return boost::math::abs( a - b );
}

template< template < typename, std::size_t > class Vec, typename T, std::size_t N >
static T meanSummedDiff( const typename std::vector< Vec< T, N > >& vecA, const typename std::vector< Vec< T, N > >& vecB )
{
	typedef typename std::vector< Vec< T, N > > container_type;
	typedef typename container_type::value_type VecType;
	
	const std::size_t n = std::distance( vecA.begin(), vecA.end() );
	const std::size_t nB = std::distance( vecB.begin(), vecB.end() );
	assert( n == nB );
	
	std::vector< VecType > vecAfter;
	vecAfter.reserve( n );
	std::transform( vecA.begin(), vecA.end(), vecB.begin(), std::back_inserter( vecAfter ), std::minus< VecType >() );
	std::vector< T > distsAfter;
	distsAfter.reserve( n );
	std::transform( vecAfter.begin(), vecAfter.end(), std::back_inserter( distsAfter ), Ubitrack::Math::Norm_2() );
	return ( std::accumulate( distsAfter.begin(), distsAfter.end(), static_cast< T > ( 0 ) ) / n );
}

template< typename T >
static T meanSummedTranslationDiff( const std::vector< Ubitrack::Math::Pose >& poseA, const std::vector< Ubitrack::Math::Pose >& poseB )
{
	const std::size_t n = std::distance( poseA.begin(), poseA.end() );
	const std::size_t nB = std::distance( poseB.begin(), poseB.end() );
	assert( n == nB );
	
	T sum = 0;
	std::vector< Ubitrack::Math::Pose >::const_iterator itB = poseB.begin();
	std::vector< Ubitrack::Math::Pose >::const_iterator itA = poseA.begin();
	std::vector< Ubitrack::Math::Pose >::const_iterator itEnd = poseA.end();
	for( ; itA<itEnd; ++itA, ++itB )
		sum += Ubitrack::Math::Norm_2() ( static_cast< Ubitrack::Math::Vector3d > ( itA->translation() - itB->translation() ) );
	return (sum/n);
}

template< typename T >
static T meanSummedAngularDiff( const std::vector< Ubitrack::Math::Pose >& poseA, const std::vector< Ubitrack::Math::Pose >& poseB )
{
	const std::size_t n = std::distance( poseA.begin(), poseA.end() );
	const std::size_t nB = std::distance( poseB.begin(), poseB.end() );
	assert( n == nB );
	
	T sum = 0;
	std::vector< Ubitrack::Math::Pose >::const_iterator itB = poseB.begin();
	std::vector< Ubitrack::Math::Pose >::const_iterator itA = poseA.begin();
	std::vector< Ubitrack::Math::Pose >::const_iterator itEnd = poseA.end();
	for( ; itA<itEnd; ++itA, ++itB )
		sum += Ubitrack::Math::Quaternion( itA->rotation() * ~(itB->rotation()) ).normalize().angle(); 
	return (sum/n);
}


#endif


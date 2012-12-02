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

#include <boost/numeric/ublas/matrix.hpp>

#include <utMath/Quaternion.h>
#include <utMath/Vector.h>
#include <math.h>

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

template< int N, typename T > 
Ubitrack::Math::Vector< N, T > randomVector( T maxVal = 100.0 )
{
	Ubitrack::Math::Vector< N, T > v;
	for ( unsigned i = 0; i < v.size(); i++ )
		v( i ) = random( T( -maxVal ), T( maxVal ) );
	return v;
}

template< class MA, class MB > 
double matrixDiff( const MA& ma, const MB& mb )
{
	double d = 0.0;
	for ( unsigned i = 0; i < ma.size1(); i++ )
		for ( unsigned j = 0; j < ma.size2(); j++ )
			d += fabs( ma( i, j ) - mb( i, j ) );
	return d / boost::numeric::ublas::norm_frobenius( ma );
}

template< class VA, class VB > 
double vectorDiff( const VA& va, const VB& vb )
{
	double d = 0.0;
	for ( unsigned i = 0; i < va.size(); i++ )
		d += fabs( va( i ) - vb( i ) );
	return d / boost::numeric::ublas::norm_2( va );
}

template< class MA, class MB > 
double homMatrixDiff( const MA& A, const MB& B )
{
	// normalize both H and Htest
	// ublas is evil. Don't try this: Htest /= ublas::norm_frobenius( Htest )
	double normA = boost::numeric::ublas::norm_frobenius( A );
	double normB = boost::numeric::ublas::norm_frobenius( B );
	
	double dp = 0.0; // sum of differences A-B
	double dm = 0.0; // sum of differences A+B (in case A ~ -B)
	for ( unsigned i = 0; i < A.size1(); i++ )
		for ( unsigned j = 0; j < A.size2(); j++ )
		{
			dp += fabs( A( i, j ) / normA - B( i, j ) / normB );
			dm += fabs( A( i, j ) / normA + B( i, j ) / normB );
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

#endif


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
 * @ingroup math
 * @file
 * Defines a function that normalizes a vector
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_VECTOR_NORMALIZE_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_VECTOR_NORMALIZE_H_INCLUDED__
 
 
 
#include <utMath/Vector.h>

namespace Ubitrack { namespace Math { namespace Optimization { namespace Function {

/**
 * A function that normalizes a vector with jacobian computation.
 * Modeled after \c UnaryFunctionPrototype
 */
class VectorNormalize
{
public:
	VectorNormalize( unsigned s )
		: m_size( s )
	{}
	
	unsigned size() const
	{ return m_size; }
	
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		typename VT2::value_type norm(  boost::numeric::ublas::norm_2( input ) );
		result = input * ( 1 / norm );
	}

	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		typename VT2::value_type norm(  boost::numeric::ublas::norm_2( input ) );
		result = input * ( 1 / norm );
		buildJacobian( input, jacobian, norm );
	}
	
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& jacobian ) const
	{
		buildJacobian( input, jacobian, boost::numeric::ublas::norm_2( input ) );
	}

protected:
	template< class VT2, class MT > 
	void buildJacobian( const VT2& input, MT& jacobian, typename VT2::value_type norm ) const
	{
		typedef typename VT2::value_type VType;
		
		// f = norm( 1 / norm( input )^3
		VType normSq = norm * norm;
		VType f( 1 / ( normSq * norm ) );

		// diagonal: ( norm( input )^2 - input( i )^2 ) * f
		for ( unsigned i = 0; i < m_size; i++ )
			jacobian( i, i ) = ( normSq - input( i ) * input( i ) ) * f;
		
		// off-diagonal: -input( i )*input( j ) * f
		for ( unsigned j = 0; j < m_size; j++ )
			for ( unsigned i = j + 1; i < m_size; i++ )
			{
				jacobian( i, j ) = - input( i ) * input( j ) * f;
				jacobian( j, i ) = jacobian( i, j );
			}
	}

	unsigned m_size;
};

}}}} // namespace Ubitrack::Math::Optimization::Function

#endif // __UBITRACK_MATH_FUNCTION_VECTOR_NORMALIZE_H_INCLUDED__
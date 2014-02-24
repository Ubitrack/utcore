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


#ifndef __UBITRACK_MATH_STOCHASTIC_MAHALANOBIS_DISTANCE_H__
#define __UBITRACK_MATH_STOCHASTIC_MAHALANOBIS_DISTANCE_H__
 
#include "Gaussian.h"

// Ubitrack
#include <utUtil/Exception.h>
#include <utUtil/StaticAssert.h>
#include "../Functors/MatrixFunctors.h"

namespace Ubitrack{ namespace Math { namespace Stochastic {


// empty struct, can be specialised further on for different types of Gaussian distributions
template< typename T >
struct MahalanobisDistance{};

template< typename T, std::size_t N >
struct MahalanobisDistance< Gaussian< T, N > >
{
	typedef T value_type;
		
protected:
	const Math::Functors::matrix_determinant m_determinant;
	const Math::Functors::matrix_inverse m_inverter;
	const Gaussian< T, N > &gaussian;
	value_type inv_covariance [ N * N ];
	bool valid;
	
public:
	MahalanobisDistance( const Gaussian< T, N >& gauss )
		: m_determinant( )
		, m_inverter( )
		, gaussian( gauss )
		, valid( false )
		{
			std::copy( gaussian.covariance, gaussian.covariance+(N*N), inv_covariance );
			
			Math::Matrix< value_type, N, N > CovMat( inv_covariance );
			const value_type det = m_determinant( CovMat );
			
			if( det != det || det == 0 ) //trick to check if the value is valid
				UBITRACK_THROW( "Cannot calculate covariance inverse, determinant is NaN." );
			
			CovMat = m_inverter( CovMat );
			std::copy( CovMat.content(), CovMat.content()+(N*N), inv_covariance );
		}
	
	
	// calculate a distance of the point
	template< typename vector_type >
	value_type operator()( const vector_type& vec ) const
	{
		// calculate a difference vector: n-diff
		value_type value[ N ];
		for( std::size_t i( 0 ); i<N; ++i )
			value[ i ] = ( vec[ i ] - gaussian.mean[ i ] );
		
		// generate a temporary result
		value_type sol_vec[ N ];
		std::fill( sol_vec, sol_vec+N, 0 );
		
		// calculate: (n-by-n covariance^-1) * n-diff
		for( std::size_t i1( 0 ); i1<N; ++i1 )
			for( std::size_t i2( 0 ); i2<N; ++i2 )
				sol_vec[ i1 ] += inv_covariance[ i1*N+i2 ] * value[ i2 ];

		// dot product of n-diff and temporary result
		value_type dot_product( 0 );
		for( std::size_t i( 0 ); i<N; ++i )
			dot_product += sol_vec[ i ] * value[ i ];
			
		return std::sqrt( dot_product );
	}
};

} } } // namespace Ubitrack::Math::Stochastic

#endif //__UBITRACK_MATH_STOCHASTIC_MAHALANOBIS_DISTANCE_H__
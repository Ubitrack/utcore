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


#ifndef __UBITRACK_EXPECTATION_MAXIMIZATION_H__
#define __UBITRACK_EXPECTATION_MAXIMIZATION_H__

// Ubitrack
#include "Gaussian.h"
#include "Weighted.h"
#include <utUtil/Exception.h>
#include <utUtil/StaticAssert.h>
#include "../Functors/MatrixFunctors.h"
 
// std
#include <vector>
#include <limits> 
#include <numeric>
#include <assert.h>
#include <algorithm>
#include <functional>

namespace Ubitrack{ namespace Math { namespace Stochastic {


/// @internal a functor template to estimate the probability of a value belonging to a gaussian distribution
template< typename GaussianType >
struct Probability
{
	typedef typename GaussianType::value_type value_type;
	static const std::size_t size = GaussianType::size;

	protected:
		const Math::Functors::matrix_determinant m_determinant;
		const Math::Functors::matrix_inverse m_inverter;
		
		const GaussianType& gaussian;
		value_type inv_covariance [ size * size ];
		value_type constant;
		
	public:
	
	Probability( const GaussianType& gaussianIn )
		: m_determinant( )
		, m_inverter( )
		, gaussian( gaussianIn )
		, constant( 0 )
		{
			std::copy( gaussian.covariance, gaussian.covariance+(size*size), inv_covariance );
			
			Math::Matrix< value_type, size, size > CovMat( inv_covariance );
			const value_type det = m_determinant( CovMat );
			
			// std::cout << "Determinant " << det << std::endl;
			
			if( det == det && std::fabs( det ) > 1e-10 ) //trick to check if the value is valid
			{
				// UBITRACK_THROW( "Cannot calculate covariance inverse, determinant is NaN or zero." );
			
				const value_type tmpconstant = 1 / std::sqrt( std::pow( static_cast< value_type >( 2.0 * 3.14159 ), size * static_cast< value_type >( 0.5 ) ) * std::fabs ( det ) );
				constant = tmpconstant;
				// if( constant != constant ) //trick to check if the value is valid
					// UBITRACK_THROW( "Cannot reliable determine constant value of probability density function, it is NaN." );
				
				CovMat = m_inverter( CovMat );
				// std::cout << CovMat << std::endl;
				std::copy( CovMat.content(), CovMat.content()+(size*size), inv_covariance );
			}
		}
	
	template< typename vector_type >
	value_type operator()( const vector_type& vec ) const
	{
		if( constant != constant || constant == 0 ) //trick to check if the value is valid
			return 0;
			// UBITRACK_THROW( "Cannot reliable determine probability of value belonging to the distribution, constant is NaN." );
			
		value_type value[ size ];
		for( std::size_t i( 0 ); i<size; ++i )
			value[ i ] = ( vec[ i ] - gaussian.mean[ i ] );
		
		value_type sol_vec[ size ];
		std::fill( sol_vec, sol_vec+size, 0 );
		
		for( std::size_t i1( 0 ); i1<size; ++i1 )
			for( std::size_t i2( 0 ); i2<size; ++i2 )
				sol_vec[ i1 ] += inv_covariance[ i1*size+i2 ] * value[ i2 ];

		value_type sum_value( 0 );
		for( std::size_t i( 0 ); i<size; ++i )
			sum_value += sol_vec[ i ] * value[ i ];
			
		const value_type return_value = constant * std::exp( -0.5 * sum_value );

		return return_value;//( return_value == return_value ) ? ( return_value > 1 ? 1 : return_value ) : 0 ;
	}
};

/// @internal a function to calculate the summarized log likelihood of many values and many probability distributions
template< typename T, typename InputIterator1, typename InputIterator2 >
T log_likelihood( const InputIterator1 ipBegin, const InputIterator1 ipEnd, const InputIterator2 ivBegin, const InputIterator2 ivLast )
{
	typedef typename std::iterator_traits< InputIterator1 >::value_type pd_type;
	typedef typename pd_type::value_type value_type;
	
	// number of values
	const std::size_t n ( std::distance( ivBegin, ivLast ) ) ;

	// initialize the final result with zeros
	std::vector< T > summen( n , 0 );
	
	// calculate for each (weighted) Gaussian
	for( InputIterator1 pdfIter( ipBegin ) ; pdfIter != ipEnd; ++pdfIter )
	{
		std::vector< T > summen_tmp;
		summen_tmp.reserve( n );
		// std::cout << "Input : " << (*pdfIter ) << "\n";
		std::transform( ivBegin, ivLast, std::back_inserter( summen_tmp ), Probability< pd_type > ( *pdfIter ) );
		std::transform( summen_tmp.begin(), summen_tmp.end(), summen_tmp.begin(), std::bind1st( std::multiplies< T >(), pdfIter->weight ) );
		std::transform( summen_tmp.begin(), summen_tmp.end(), summen.begin(), summen.begin(), std::plus< T >() );
	}
	std::transform( summen.begin(), summen.end(), summen.begin(), std::ptr_fun< T, T >( std::log ) );
	const T loglikelihood = std::accumulate( summen.begin(), summen.end(), static_cast< T >( 0 ) );
	return loglikelihood /  static_cast< T >( n );
}


/**
 * Expectation Maximization
 *
 * This is a template-framework to carry out the expectation maximization algorithm.
 * 
 * @tparam InputIterator defines the type of input values
 * @tparam OutputIterator defines the type of output values (Probability Distribution)
 * @param itBegin
 * @param itEnd
 * @param itBeginGauss
 * @param itEndGauss
 */

template< typename InputIterator, typename OutputIterator >
typename std::iterator_traits< OutputIterator >::value_type::value_type expectation_maximization( 
	const InputIterator itBegin, const InputIterator itEnd
	, OutputIterator itBeginGauss, OutputIterator itEndGauss )
{
	/// @todo need to check if this type has an susbscript operator operator[]
	// the type of the multivariate data, e.g. float/double/int/char...
	typedef typename std::iterator_traits< InputIterator >::value_type value_type;
	
	//type of the probability density function
	typedef typename std::iterator_traits< OutputIterator >::value_type pdf_type;
	
	//some general abbreviations, needed in the algorithm
	typedef typename pdf_type::value_type numeric_type;
	typedef typename std::vector< numeric_type > numeric_container_type;
	typedef typename numeric_container_type::iterator NumIterator;
	
	//determine the number of elements/ pdfs
	const std::size_t n = std::distance( itBegin, itEnd );
	const std::size_t k_cluster ( std::distance( itBeginGauss, itEndGauss ) );
	
	// you never should have less values than expected clusters:
	assert( n > k_cluster );
		
	const numeric_type threshold( static_cast< numeric_type > ( 1e-5 ) );
	const std::size_t max_iter( 100 );
	
	//assign the norms
	// norm : norm_{1} ,.., norm_{n}
	std::vector< numeric_type > norms( n , 0 );
	// norms.reserve( n );
	// norms.assign( n, 0 );
	
	
	//assign the gamma values
	// gamma : g_{1,1} ,.., g_{n,1} ,..., g_{1,k} ,..., g_{n,k}
	std::vector< numeric_type > gammas( n*k_cluster, 0 );
	// gammas.reserve( n*k_cluster );
	// gammas.assign( n*k_cluster, 0 );
		
	numeric_type likelihood = log_likelihood< numeric_type >( itBeginGauss, itEndGauss, itBegin, itEnd );
	// std::cout << "1st Likelihood " << likelihood << std::endl;	
	
	std::size_t i( 0 );
	for( ; i< max_iter; ++i )
	{
		
		
		NumIterator gammaIter( gammas.begin() );
		
		// Expectation Step:
		for( OutputIterator pdfIter( itBeginGauss ); pdfIter< itEndGauss; ++pdfIter )
		{
			const Probability< pdf_type > pdf( *pdfIter );
			NumIterator normIter( norms.begin() );
			for( InputIterator valueIter( itBegin ); valueIter != itEnd; ++valueIter, ++normIter, ++gammaIter  )
			{
				(*normIter) += (*gammaIter) = pdfIter->weight * pdf( *valueIter );
				// std::cout << "Probability of " << *valueIter << " at " << *pdfIter << ":\n" << pdf( *valueIter ) << "\n";
			}
		}
			
		
		// const numeric_type summe_all = std::accumulate( gammas.begin(), gammas.end(), static_cast< numeric_type >( 0 ) );
		//std::transform( norms.begin(), norms.end(), norms.begin(), std::bind1st( std::divides< numeric_type >(), 1./static_cast< numeric_type > ( n ) ) );
		
		gammaIter = gammas.begin();
		for( OutputIterator pdfIter( itBeginGauss ); pdfIter != itEndGauss ; std::advance( gammaIter, n ), ++pdfIter )
		{
			// std::transform( gammaIter, gammaIter+n, norms.begin(), gammaIter, std::multiplies< numeric_type >() );
			std::transform( gammaIter, gammaIter+n, norms.begin(), gammaIter, std::divides< numeric_type >() );
			const numeric_type summe = std::accumulate( gammaIter, gammaIter+n, static_cast< numeric_type >( 0 ) );
			if( summe == 0 )
			{
				pdfIter->weight = 0;
				// std::cout << "Sum of gammas resulted in zero " << std::endl;
				continue;
			}
			

			// std::cout << " Distance " << std::distance( gammas.begin(), gammaIter ) << "\n";
			std::transform( gammaIter, gammaIter+n, gammaIter, std::bind2nd( std::multiplies< numeric_type >(), 1./summe ) );
			
			// std::cout << " Sum for the PDFS " << std::accumulate( gammaIter, gammaIter+n, static_cast< numeric_type >( 0 ) )  << std::endl;
			pdfIter->weight = summe / n ;
			if( pdfIter->weight != pdfIter->weight )
				pdfIter->weight = 0;
			else
				estimate_gaussian( itBegin, itEnd, gammaIter, *pdfIter );
			// std::cout << "Gaussian " <<  (*pdfIter) << std::endl;
		}
		// {
			// numeric_type sum = 0;
			// for( OutputIterator pdfIter( itBeginGauss ); pdfIter != itEndGauss ; ++pdfIter )
				// sum += pdfIter->weight;
			// std::cout << "Sum of all weights (should correspond to one ): " << sum <<  std::endl;
		// }
		
		//check convergence criteria
		const numeric_type newLikelihood = log_likelihood< numeric_type >( itBeginGauss, itEndGauss, itBegin, itEnd );
		if ( std::abs( likelihood - newLikelihood ) <  ( threshold * std::fabs( likelihood ) ) )
			return newLikelihood;
			
		if( newLikelihood != newLikelihood )
			return likelihood;
		likelihood = newLikelihood;
		// std::cout << "Likelihood " << likelihood << std::endl;
		
		std::fill( norms.begin(), norms.end(), 0 );
	}
	
	// std::cout << "Gaussians after " << i << " iterations\n";
	// for( OutputIterator pdfIter( itBeginGauss ); pdfIter != itEndGauss ; ++pdfIter )
		// std::cout << (pdfIter->gaussian) << std::endl;
	return likelihood;
};

} } } // namespace Ubitrack::Math::Stochastic

#endif //__UBITRACK_EXPECTATION_MAXIMIZATION_H__
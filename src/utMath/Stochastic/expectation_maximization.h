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
 
// std
#include <vector>
#include <limits> 
#include <numeric>
#include <assert.h>
#include <algorithm>
#include <functional>

// Ubitrack
#include <utUtil/Exception.h>
#include <utUtil/StaticAssert.h>
#include "../Functors/MatrixFunctors.h"

namespace Ubitrack{ namespace Math { namespace Stochastic {

/// @internal a generic type that can represented a weighted anything -> own header?
template< typename any_class, typename T >
struct Weighted
	: public any_class
{

	typedef T								weight_type;
	typedef any_class						base_type;
	typedef typename any_class::value_type	value_type;
	
	T weight;
		
	Weighted( )
		: base_type( )
		, weight ( 0 )
	{};
	
	
	bool operator< ( const Weighted< any_class, T >& other ) const 
	{
		return (this->weight < other.weight);
	}
		
};

/// @internal a generic Gaussian distribution type -> own header
template< typename T, std::size_t N >
struct Gaussian
{
	/** simple typedef to type of Gaussian structure. */
	typedef T value_type;
	typedef std::size_t size_type;
	static const size_type size = N; 
		
	/** the mean value */
	T mean[ N ];
	
	/** the covariance matrix, upper triangle equals the lower triangle */
	T covariance [ N * N ];
	
	/** the root of the sum of squared diagonal entries of the covariance */
	T variance;
	
	/** the sum of squared diagonal entries of the covariance */
	T variance2;
};

/** @internal overrides the stream output to have nicely aligned data */
template< typename T, std::size_t N >
std::ostream& operator<<( std::ostream& s, const Gaussian< T, N >& gauss )
{
	for( std::size_t i1( 0 ); i1<N; ++i1 )
	{
		s << std::setfill(' ')
		<< std::setw(10)
		<< std::fixed
		<< std::setprecision(2)
		<< gauss.mean[ i1 ] << " [ ";
		for( std::size_t i2( 0 ); i2<N; ++i2 )
			s << std::setw( 10 )
			<< std::fixed
			<< std::setprecision(4)
			<< gauss.covariance[ i1*N+i2 ] ;
		s << " ]" << std::endl;
	}
	return s;
}

/// @internal a function that estimates a Gaussian distribution (new header?)
template< typename T, std::size_t N, typename InputIterator1, typename InputIterator2 >
bool estimate_gaussian( const InputIterator1 itBegin, const InputIterator1 itEnd, const InputIterator2 itWeights, Gaussian< T, N > &gaussian )
{
	const std::size_t n = std::distance( itBegin, itEnd );
	
	if( n == 0 )// no value exit early
		return false;
		
	// set the covariance to zero elements
	std::fill( gaussian.covariance, gaussian.covariance+(N*N), static_cast< T >( 0 ) );
	if( n == 1 ) //only one value, simple case
	{
		for( std::size_t i( 0 ); i < N; ++i )
			gaussian.mean[ i ] = (*itBegin)[ i ];
		gaussian.variance = gaussian.variance2 = 0;
		return true;
	}
	// since we got that far also zero out the mean value
	std::fill( gaussian.mean, gaussian.mean+N, static_cast< T >( 0 ) );
	// Calculate (weighted) meanvalue
	{
		InputIterator2 weightIter( itWeights );
		for( InputIterator1 valIter ( itBegin ); valIter != itEnd; ++valIter, ++weightIter )
			for( std::size_t i( 0 ); i < N; ++i )
				gaussian.mean[ i ] += (*weightIter) * (*valIter)[ i ];
	}
		
	// Calculate (weighted) covariance
	{
		InputIterator2 weightIter( itWeights );
		for( InputIterator1 valIter ( itBegin ); valIter != itEnd; ++valIter, ++weightIter )
			for( std::size_t i1( 0 ); i1 < N; ++i1 )
				for( std::size_t i2( i1 ); i2 < N; ++i2 )
					gaussian.covariance[ i2*N+i1 ] = gaussian.covariance[ i1*N+i2 ] += ( *weightIter ) * ((*valIter)[ i1 ] - gaussian.mean[ i1 ] ) * ( (*valIter)[ i2 ] - gaussian.mean[ i2 ] );
	}
	
	{
		gaussian.variance2 = 0;
		//sum up the squared diagonal entries
		for( std::size_t i( 0 ); i < N; ++i )	
			gaussian.variance2 += ( gaussian.covariance[i*N+i] * gaussian.covariance[i*N+i] );
			
		gaussian.variance = std::sqrt( gaussian.variance2 );
	}
	return true;
};

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
		{
			std::copy( gaussian.covariance, gaussian.covariance+(size*size), inv_covariance );
			
			Math::Matrix< value_type, size, size > CovMat( inv_covariance );
			const value_type det = m_determinant( CovMat );
			
			if( det != det ) //trick to check if the value is valid
				UBITRACK_THROW( "Cannot calculate covariance inverse, determinant is NaN." );
				
			// const value_type sqrtDet = std::sqrt( std::fabs( det ) );
			
			const value_type tmpconstant = 1.0 / std::sqrt( std::pow( 2.0 * 3.14159, static_cast< value_type >( size ) ) * std::fabs( det ) );
			
			constant = tmpconstant;
			if( constant != constant ) //trick to check if the value is valid
				UBITRACK_THROW( "Cannot reliable determine constant value of probability density function, it is NaN." );
			
			CovMat = m_inverter( CovMat );
			// std::cout << CovMat << std::endl;
			std::copy( CovMat.content(), CovMat.content()+(size*size), inv_covariance );
		}
	
	template< typename vector_type >
	value_type operator()( const vector_type& vec ) const
	{
	
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
			
		const value_type return_value = constant * std::exp(  - 0.5 * sum_value );

		return ( return_value == return_value ) ? ( return_value > 1 ? 1 : return_value ) : 0 ;
	}
};

/// @internal a function to calculate the summarized log likelihood of many values and many probability distributions
template< typename T, typename InputIterator1, typename InputIterator2 >
T log_likelihood( const InputIterator1 ipBegin, const InputIterator1 ipEnd, const InputIterator2 ivBegin, const InputIterator2 ivLast )
{
	typedef typename std::iterator_traits< InputIterator1 >::value_type pd_type;
	typedef typename pd_type::value_type value_type;
	
	
	const std::size_t n ( std::distance( ivBegin, ivLast ) ) ;

	std::vector< T > summen;
	summen.assign( n , 0 );
	for( InputIterator1 pdfIter( ipBegin ) ; pdfIter != ipEnd; ++pdfIter )
	{
		std::vector< T > summen_tmp;
		summen_tmp.reserve( n );
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
	std::vector< numeric_type > norms;
	norms.reserve( n );
	norms.assign( n, 0 );
	
	//assign the gamma values
	// gamma : g_{1,1} ,.., g_{n,1} ,..., g_{1,k} ,..., g_{n,k}
	std::vector< numeric_type > gammas;
	gammas.reserve( n*k_cluster );
	gammas.assign( n*k_cluster, 0 );
		
	numeric_type likelihood = log_likelihood< numeric_type >( itBeginGauss, itEndGauss, itBegin, itEnd );
// std::cout << "1st Likelihood " << likelihood << std::endl;	
	
	std::size_t i( 0 );
	for( ; i< max_iter; ++i )
	{
		
		std::fill( norms.begin(), norms.end(), 0 );
		NumIterator gammaIter( gammas.begin() );
		
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
		std::transform( norms.begin(), norms.end(), norms.begin(), std::bind1st( std::divides< numeric_type >(), 1./static_cast< numeric_type > ( n ) ) );
		
		gammaIter = gammas.begin();
		for( OutputIterator pdfIter( itBeginGauss ); pdfIter != itEndGauss ; std::advance( gammaIter, n ), ++pdfIter )
		{
			std::transform( gammaIter, gammaIter+n, norms.begin(), gammaIter, std::multiplies< numeric_type >() );
			const numeric_type summe = std::accumulate( gammaIter, gammaIter+n, static_cast< numeric_type >( 0 ) );
			if( summe == 0 )
			{
				pdfIter->weight = 0;
				std::cout << " Gamma resulted in zero " << std::endl;
				continue;
			}
			
			// std::cout << " Sum for the PDFS " << summe << std::endl;
			
			std::transform( gammaIter, gammaIter+n, gammaIter, std::bind2nd( std::multiplies< numeric_type >(), 1./summe ) );
			estimate_gaussian( itBegin, itEnd, gammaIter, *pdfIter );
			
			pdfIter->weight = summe ;
			// std::cout << " Summe gewicht  " << summe <<  std::endl;
			// pdfIter->weight = summe / summe_all;
		}
		
		
		//check convergence criteria
		const numeric_type newLikelihood = log_likelihood< numeric_type >( itBeginGauss, itEndGauss, itBegin, itEnd );
		if ( std::abs( likelihood - newLikelihood ) <  ( threshold * std::fabs( likelihood ) ) )
			return newLikelihood;
			
		likelihood = newLikelihood;
		// std::cout << "Likelihood " << likelihood << std::endl;
	}
	
	// std::cout << "Gaussians after " << i << " iterations\n";
	// for( OutputIterator pdfIter( itBeginGauss ); pdfIter != itEndGauss ; ++pdfIter )
		// std::cout << (pdfIter->gaussian) << std::endl;
	return likelihood;
};

} } } // namespace Ubitrack::Math::Stochastic

#endif //__UBITRACK_EXPECTATION_MAXIMIZATION_H__
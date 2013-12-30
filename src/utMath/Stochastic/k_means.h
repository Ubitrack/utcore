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
 * @ingroup math stochasic
 * @file
 *
 * Several functions necessary for k-means calculations including the k-means algorithm.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_K_MEANS_H__
#define __UBITRACK_K_MEANS_H__
 

// Ubitrack
// #include <../Vector.h>
#include <utUtil/StaticAssert.h>
#include "../Functors/VectorFunctors.h"
#include "../Random/Scalar.h" // Randomness needed for kmeans++
#include "../Geometry/container_traits.h"


// std
#include <vector>
#include <limits> 
#include <numeric>
#include <algorithm>
#include <functional>

namespace {

// small template to support k-means
template< typename InputIterator  >
struct assign_indices
{
	// typedef typename std::iterator_traits< ForwardIterator >::value_type vector_type;
	const InputIterator iBegin;
	const InputIterator iEnd;
	
	assign_indices( const InputIterator first , const InputIterator last )
		: iBegin( first )
		, iEnd( last )
		{}
	

	template< template< typename, std::size_t > class VecType, typename T, std::size_t N >
	std::size_t operator()( const typename VecType< T, N >& vec ) const
	{
		std::size_t k( 0 );
		InputIterator iter = iBegin;
		T d = Ubitrack::Math::Functors::Norm_1< T, N >()( vec - (*iter++) );
		
		for( std::size_t k_now( 1 ) ; iter<iEnd; ++iter, ++k_now )
		{
			const T d_new = Ubitrack::Math::Functors::Norm_1< T, N >()( vec - (*iter) );
			if( d_new < d )
			{
				d = d_new;
				k = k_now;
			}
		}
		return k;
	}
};

// small template to support k-means
template< typename InputIterator >
struct assign_to_mean
{
	InputIterator iBegin;
	
	assign_to_mean( InputIterator first )
		: iBegin( first ) {}

	template< typename VecType, typename iType >
	void operator()( const VecType& vec, const iType index )
	{
		*(iBegin + index) += vec;
	}
};

// a helper function
template< typename T1, typename T2, typename BinaryFunction >
BinaryFunction k_means_accumulate( const T1 pFirst1, const T1 pLast, const T2 pFirst2, BinaryFunction op )
{
	T1 first1 = pFirst1;
	T2 first2 = pFirst2;
	for( ; first1 != pLast; ++first1, ++first2 )
		op( *first1, *first2 );

	return op;
};

} // namespace anonymous


namespace Ubitrack{ namespace Math { namespace Stochastic {
 
/**
 * @ingroup math stochastic
 * @brief picks the first k elements of a given range.
 *
 * This function can be applied to nearly any container structure
 * providing access to the single container elements via forward iterators.
 * Although designed for stl-containers * (e.g. \c std::vector,
 * \c std::list or \c std::set, etc ) the algorithm is not limited
 * to these ones but can be used via any pointer to datastructures.
 *
 * Example use case:\n
 * std::vector< Vector3d > points3d; // <- class to pick some \c k elements \n
 * std::vector< Vector3d > points3dOut; // <- will be filled with values, storage can be allocated with \c reserve() \n
 * copy_greedy( points3d.begin(), points3d.end(), k, std::back_inserter( points3dOut ) );\n
 * 
 * @tparam ForwardIterator1 type of forward iterator to container of input elements
 * @tparam ForwardIterator2 type of the output iterator to container for selected input elements
 * @param iBegin \c iterator pointing to first element in the input container/storage class of the elements ( usually \c begin() )
 * @param iEnd \c iterator pointing behind the last element in the input container/storage class of the elements ( usually \c end() )
 * @param n_cluster a value that signs how many elements should be selected and put to the output iterator
 * @param itSelected output \c iterator pointing to first element in container/storage class for storing the selected elements ( usually \c begin() or \c std::back_inserter(container) )
 */ 
template< typename ForwardIterator, typename OutputIterator >
void copy_greedy( ForwardIterator iBegin, ForwardIterator iEnd, const std::size_t n_cluster, OutputIterator itSelected )
{
	typedef typename std::iterator_traits< ForwardIterator >::value_type vector_type;
	for( std::size_t k( 0 ); k < n_cluster; ++k )
		(*itSelected++) = ( *iBegin++ );
};



/// internal function only, please see other k-means for explanation
template<  typename ForwardIterator, typename OutputIterator, typename IndicesIterator, template < typename, std::size_t > class VecType, typename T, std::size_t N >
T k_means( const ForwardIterator iBegin, const ForwardIterator iEnd, const std::size_t n_cluster, OutputIterator itOut, IndicesIterator indicesOut, typename VecType< T, N > )
{
	// some typedefs regarding the input vector type
	typedef typename std::iterator_traits< ForwardIterator >::value_type vector_type;
	typedef typename vector_type::value_type value_type;
	
	// some typedefs regarding the indices container
	typedef typename Ubitrack::Util::container_traits< IndicesIterator >::value_type indices_type;
	typedef typename std::vector< indices_type > indices_container_type;
	
	typedef typename std::vector< vector_type > mean_container_type;
	typedef typename mean_container_type::iterator mean_type_iterator;
	
	// we use a quadratic type for epsilon, such that we do not need to calculate square roots later
	const value_type epsilon = std::pow( static_cast< value_type > (1e-02), 2 );
	const std::size_t max_iter = 100; // names say all, I hope
	
	// std::cout << "Vector Type " << typeid( vector_type ).name() << " of dimension=" << N << " and type=" << typeid( T ).name() << "\n";	
	
	const std::size_t n = std::distance( iBegin, iEnd );
	assert( n > n_cluster );
	// std::cout << "k-means called with " << n << " elements to determine " << n_cluster << " clusters\n";
	
	// the sum of the different means, saved separately for each dimension.
	mean_container_type means;
	means.reserve( n_cluster );
	copy_greedy( iBegin, iEnd, n_cluster, std::back_inserter( means ) );

	//assign indices for the first time
	indices_container_type indices;
	indices.reserve( n_cluster );
	std::transform( iBegin, iEnd, std::back_inserter( indices ), assign_indices< mean_type_iterator >( means.begin(), means.end() ) );
	
	
	std::size_t i( 0 );
	value_type diff_error;
	for( ; i<max_iter; ++i )
	{
		// set some "fresh" (empty) mean values 
		mean_container_type means_temp;
		means_temp.reserve( n_cluster );
		means_temp.assign( n_cluster, vector_type::zeros() );
				
		//accumulate the means from the clusters 
		k_means_accumulate( iBegin, iEnd, indices.begin(), assign_to_mean< mean_type_iterator >( means_temp.begin() ) );
		
		// reset the amount of gathered values, could be done nicer...
		for( std::size_t k = 0; k<n_cluster; ++k )		
			means_temp[ k ] /= std::count ( indices.begin(), indices.end(), k );
		
		// mean_container_type means_diff;
		// means_diff.reserve( n_cluster );
		// std::transform( means_temp.begin(), means_temp.end(), means.begin(), std::back_inserter( means_diff ), std::minus< vector_type >() );
		// generate the difference vectors
		std::transform( means_temp.begin(), means_temp.end(), means.begin(), means.begin(), std::minus< vector_type >() );
		
		// calculate the summarized difference
		std::vector< value_type > norms1;
		norms1.reserve( n_cluster );
		std::transform( means.begin(), means.end(),  std::back_inserter( norms1 ), Ubitrack::Math::Functors::Norm_1< T, N >() );
		diff_error = std::accumulate( norms1.begin(), norms1.end(), static_cast< value_type >( 0 ) ) / n_cluster;
		
		// reset the means values
		means.assign( means_temp.begin(), means_temp.end() );

		// finally assign the indices to the corresponding clusters for the new loop
		std::transform( iBegin, iEnd, indices.begin(), assign_indices< mean_type_iterator >( means.begin(), means.end() ) );
		
		if( diff_error < epsilon )
			break;
	}
	
	
	// copy the resulting mean values to output iterator
	std::copy( means.begin(), means.end(), itOut );
	// and copy the indices to the corresponding output iterator
	std::copy( indices.begin(), indices.end(), indicesOut );
	
	return diff_error;
};

// /// overloaded function, see other k-means for detailed description
// template< typename InputIterator, typename OutputIterator >
// void k_means( const InputIterator iBeginValues, const InputIterator iEndValues, const std::size_t n_cluster, OutputIterator itIndices )
// {
	// typedef typename std::iterator_traits< ForwardIterator >::value_type VecType;
	// std::vector< VecType > means;
	// k_means( iBeginValues, iEndValues, n_cluster, std::back_inserter( means ), itIndices, *iBeginValues );
// };

/**
 * @ingroup math stochastic
 * @brief determines k centroids of clusters from a set of elements
 *
 * This function can be applied to nearly any container structure
 * providing access to the single container elements via forward iterators.
 * Although designed for stl-containers * (e.g. \c std::vector,
 * \c std::list or \c std::set, etc ) the algorithm is not limited
 * to these ones but can be used via any pointer to data-structures.
 *
 * Example use case:\n
 * std::vector< Vector3d > points3d; // <- class to pick some \c k elements \n
 * std::size_t k; // <- defines the number of clusters to be determined \n
 * std::vector< Vector3d > points3dOut; // <- will be filled with the cluster centroids, storage can be allocated with \c reserve() \n
 * std::vector< std::size_t > indices; // <- will sign to with cluster a point belongs corresponds to the elements from the input iterators, storage can be allocated with \c reserve() \n
 * k_means( points3d.begin(), points3d.end(), k, std::back_inserter( points3dOut ), std::back_inserter( indices ) );\n
 * 
 * @tparam ForwardIterator type of forward iterator to container of input elements
 * @tparam OutputIterator1 type of the output iterator to container of the clusters
 * @tparam OutputIterator2 type of the output iterator to container of the indices
 * @param iBeginValues \c iterator pointing to first element in the input container/storage class of the input elements ( usually \c begin() )
 * @param iEndValues \c iterator pointing behind the last element in the input container/storage class of the input elements ( usually \c end() )
 * @param n_cluster a value that signs the amount of clusters and corresponding centroids that are expected within the input elements
 * @param iBeginMean output \c iterator pointing to first element in container/storage class for storing the determined centroids ( usually \c begin() or \c std::back_inserter(container) )
 * @param itIndices output \c iterator pointing to first element in container/storage class for storing the indices to the centroids of the corresponding input elements ( usually \c begin() or \c std::back_inserter(container) )
 */ 
template< typename InputIterator, typename OutputIterator1, typename OutputIterator2 >
void k_means( const InputIterator iBeginValues, const InputIterator iEndValues, const std::size_t n_cluster, OutputIterator1 itCentroids, OutputIterator2 itIndices )
{
	k_means( iBeginValues, iEndValues, n_cluster, itCentroids, itIndices, *iBeginValues );
};

} } } // namespace Ubitrack::Math::Stochastic

#endif //__UBITRACK_K_MEANS_H__
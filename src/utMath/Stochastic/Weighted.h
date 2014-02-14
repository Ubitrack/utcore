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
 * @ingroup math stochastic
 * @file
 *
 * A template class to realize a weighted version of (nearly) anything.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_MATH_WEIGHTED_H__
#define __UBITRACK_MATH_WEIGHTED_H__
 

// std
#include <iosfwd>


namespace Ubitrack{ namespace Math { namespace Stochastic {

/**
 * Weighted template class
 *
 * This template class can represent any type of weighted class.
 * It is helpful to represent Gaussian Mixture Models, Particle Filters
 * and similar probabilistic frameworks.
 *
 * This class provides some overloaded \c boolean \c comparison \c operators 
 * ( "<", ">", "==" and "!=" ) such that stl-algorithms ( like \c std::sort )
 * can be easily applied to stl-containers of weighted types.
 *
 * @tparam classType the type of the (weighted) class
 * @tparam weightType built-in type of the weight ( e.g. \c double or \c float ), defaults to value_type of the classType
 */
template< typename classType, typename weightType = typename classType::value_type >
struct Weighted
	: public classType
{

	typedef weightType							weight_type;
	typedef classType							base_type;
	typedef Weighted< base_type, weight_type >	self_type;
	typedef typename classType::value_type		value_type;
	
	weight_type weight;
		
	Weighted( )
		: base_type( )
		, weight ( 0 )
	{};
	
	Weighted( const base_type& base )
		: base_type( base )
		, weight ( 0 )
	{};
	
	Weighted( const base_type& base, const weight_type w )
		: base_type( base )
		, weight ( w )
	{};
	
	
	bool operator< ( const self_type& other ) const 
	{
		return (this->weight < other.weight);
	}
	
	bool operator> ( const self_type& other ) const 
	{
		return (this->weight > other.weight);
	}
	
	bool operator== ( const self_type& other ) const 
	{
		return (this->weight == other.weight);
	}
	
	bool operator!= ( const self_type& other ) const 
	{
		return (this->weight != other.weight);
	}
		
};

} } } // namespace Ubitrack::Math::Stochastic

#endif //__UBITRACK_MATH_WEIGHTED_H__


/** stream output operator for any weighted type */
// template< typename classType, typename weightType >
// std::ostream& operator<<( std::ostream& s, const Ubitrack::Math::Stochastic::Weighted< classType, weightType >& weighted )
// {
	// s << "weight=" << weighted.weight << "\n";
	// s << static_cast< classType >( weighted );
	// return s;
// }
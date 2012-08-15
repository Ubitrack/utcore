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
 * @file
 * Function that performs a linear time update on a vector.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 


#ifndef __UBITRACK_TRACKING_FUNCTION_LINEARTIMEUPDATE_H_INCLUDED__
#define __UBITRACK_TRACKING_FUNCTION_LINEARTIMEUPDATE_H_INCLUDED__

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace Ubitrack { namespace Tracking { namespace Function {

/**
 * Updates a vector assuming constant velocity, acceleration or higher order derivatives
 */
class LinearTimeUpdate
{
public:
	/**
	 * constructor.
	 * @param deltaTime time to forward in seconds
	 * @param vectorSize size of the vector to be updated
	 * @param order number of derivatives to take into account
	 */
	LinearTimeUpdate( double deltaTime, unsigned vectorSize, unsigned order = 1 )
		: m_deltaTime( deltaTime )
		, m_size( vectorSize )
		, m_order( order )
	{}

	unsigned size() const
	{ return m_size; }
	
	/**
	 * @param result vector to put the result in
	 * @param input vector with (vectorSize * (1 + \c order)) elements containing the input vector and \c order derivatives
	 * @param jacobian size-by-(size*(1+order)) matrix to put the jacobian into
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		namespace ublas = boost::numeric::ublas;

		result = ublas::subrange( input, 0, m_size );
		ublas::subrange( jacobian, 0, m_size, 0, m_size ) = 
			ublas::identity_matrix< typename MT::value_type >( m_size, m_size );
		
		typename MT::value_type t = 1;
		for ( unsigned i = 1; i <= m_order; i++ )
		{
			t *= static_cast< typename MT::value_type >( m_deltaTime ) / i;
			
			result += t * ublas::subrange( input, m_size * i, m_size * (i+1) );
			ublas::subrange( jacobian, 0, m_size, m_size * i, m_size * (i+1) ) = 
				t * ublas::identity_matrix< typename MT::value_type >( m_size, m_size );
		}
	}

protected:
		double m_deltaTime;
		unsigned m_size;
		unsigned m_order;
};

} } } // namespace Ubitrack::Tracking::Function

#endif

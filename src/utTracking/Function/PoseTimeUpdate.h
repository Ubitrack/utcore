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
 * Function that performs a time update on a pose and its derivatives.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 


#ifndef __UBITRACK_TRACKING_FUNCTION_POSETIMEUPDATE_H_INCLUDED__
#define __UBITRACK_TRACKING_FUNCTION_POSETIMEUPDATE_H_INCLUDED__

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "LinearTimeUpdate.h"
#include "QuaternionTimeUpdate.h"

namespace Ubitrack { namespace Tracking { namespace Function {

/**
 * Updates a pose assuming constant velocity, acceleration or higher order derivatives in both 
 * position and orientation
 */
class PoseTimeUpdate
{
public:
	/**
	 * constructor.
	 * @param deltaTime time to forward in seconds
	 * @param vectorSize size of the vector to be updated
	 * @param posOrder number of position derivatives to take into account. Choose -1 for no position.
	 * @param rotOrder number of rotation derivatives to take into account. Choose -1 for no rotation
	 */
	PoseTimeUpdate( double deltaTime, int posOrder, int rotOrder )
		: m_deltaTime( deltaTime )
		, m_posOrder( posOrder )
		, m_rotOrder( rotOrder )
	{}

	unsigned size() const
	{ return 6 + 3 * m_posOrder + 3 * m_rotOrder + ( m_rotOrder >= 0 ? 1 : 0 ); }
	
	/**
	 * Updates the vector and computes the jacobian.
	 * Both input and output vectors must be of size 7 + 3 * (rotOrder+posOrder). The layout must be
	 *  p, p', p'', ..., q, q', q'', ...
	 * @param result vector to put the result in
	 * @param input input vector 
	 * @param jacobian matrix to put the jacobian into
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		namespace ublas = boost::numeric::ublas;
		const int rs = 3 + 3 * m_posOrder; // start of rotation
		const int ts = size(); // total size
		
		// update position and its derivatives
		for ( int i = 0; i <= m_posOrder; i++ )
		{
			const int s = 3 * i;
			ublas::vector_range< VT1 > resultSubRange( result, ublas::range( s, s + 3 ) );
			ublas::matrix_range< MT > jacobianSubRange( jacobian, ublas::range( s, s + 3 ), ublas::range( s, rs ) );
			LinearTimeUpdate( m_deltaTime, 3, m_posOrder - i )
				.evaluateWithJacobian( resultSubRange, ublas::subrange( input, s, rs ), jacobianSubRange );
			ublas::subrange( jacobian, s, s + 3, 0, s ) = 
				ublas::zero_matrix< typename MT::value_type >( 3, s );
			ublas::subrange( jacobian, s, s + 3, rs, ts ) = 
				ublas::zero_matrix< typename MT::value_type >( 3, ts - rs );
		}
		
		// update the rotation quaternion
		if ( m_rotOrder >= 0 )
		{
			ublas::vector_range< VT1 > resultSubRange( result, ublas::range( rs, rs + 4 ) );
			ublas::matrix_range< MT > jacobianSubRange( jacobian, ublas::range( rs, rs + 4 ), ublas::range( rs, ts ) );
			QuaternionTimeUpdate( m_deltaTime, m_rotOrder )
				.evaluateWithJacobian( resultSubRange, ublas::subrange( input, rs, ts ), jacobianSubRange );
			ublas::subrange( jacobian, rs, rs + 4, 0, rs ) = 
				ublas::zero_matrix< typename MT::value_type >( 4, rs );
		}
		
		// update the quaternion derivatives
		for ( int i = 1; i <= m_rotOrder; i++ )
		{
			const int s = rs + 1 + 3 * i;
			ublas::vector_range< VT1 > resultSubRange( result, ublas::range( s, s + 3 ) );
			ublas::matrix_range< MT > jacobianSubRange( jacobian, ublas::range( s, s + 3 ), ublas::range( s, ts ) );
			LinearTimeUpdate( m_deltaTime, 3, m_rotOrder - i )
				.evaluateWithJacobian( resultSubRange, ublas::subrange( input, s, ts ), jacobianSubRange );
			ublas::subrange( jacobian, s, s + 3, 0, s ) = 
				ublas::zero_matrix< typename MT::value_type >( 3, s );
		}

	}

protected:
	double m_deltaTime;
	int m_posOrder;
	int m_rotOrder;
};

} } } // namespace Ubitrack::Tracking::Function

#endif

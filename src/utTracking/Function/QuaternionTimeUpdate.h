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
 * Function that performs a time update on a quaternion.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 


#ifndef __UBITRACK_TRACKING_FUNCTION_QUATERNIONTIMEUPDATE_H_INCLUDED__
#define __UBITRACK_TRACKING_FUNCTION_QUATERNIONTIMEUPDATE_H_INCLUDED__

#include <utMath/Quaternion.h>
#include <utMath/RotationVelocity.h>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace Ubitrack { namespace Tracking { namespace Function {

/**
 * Updates a quaternion assuming constant angular velocity
 */
class QuaternionTimeUpdate
{
public:
	/**
	 * constructor.
	 * @param deltaTime time to forward in seconds
	 * @param order numer of derivatives to take into account
	 * TODO: only order=1 works correctly
	 */
	QuaternionTimeUpdate( double deltaTime, unsigned order = 1 )
		: m_deltaTime( deltaTime )
		, m_order( order )
	{}

	unsigned size() const
	{ return 4; }
	
	/**
	 * @param result 4-vector to put the result quaternion in
	 * @param input vector with (4 + 3 * \c order) elements containing the input quaternions and \c order derivatives
	 * @param jacobian 4-by-(4+3*order) matrix to put the jacobian into
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		namespace ublas = boost::numeric::ublas;
		
		Math::Quaternion q( Math::Quaternion::fromVector( input ) );
		Math::Quaternion r( q );
		double t = m_deltaTime;
		for ( unsigned i = 0; i < m_order; i++ )
		{
			Math::RotationVelocity v( ublas::subrange( input, 4 + 3 * i, 4 + 3 * ( i + 1 ) ) );
			r = r * v.integrate( t );
			t *= m_deltaTime;
			
			jacobian( 0, 4 + i * 3 + 0 ) =  t * q.w() / 2;
			jacobian( 0, 4 + i * 3 + 1 ) = -t * q.z() / 2;
			jacobian( 0, 4 + i * 3 + 2 ) =  t * q.y() / 2;
			
			jacobian( 1, 4 + i * 3 + 0 ) =  t * q.z() / 2;
			jacobian( 1, 4 + i * 3 + 1 ) =  t * q.w() / 2;
			jacobian( 1, 4 + i * 3 + 2 ) = -t * q.x() / 2;

			jacobian( 2, 4 + i * 3 + 0 ) = -t * q.y() / 2;
			jacobian( 2, 4 + i * 3 + 1 ) =  t * q.x() / 2;
			jacobian( 2, 4 + i * 3 + 2 ) =  t * q.w() / 2;
			
			jacobian( 3, 4 + i * 3 + 0 ) = -t * q.x() / 2;
			jacobian( 3, 4 + i * 3 + 1 ) = -t * q.y() / 2;
			jacobian( 3, 4 + i * 3 + 2 ) = -t * q.z() / 2;
		}
		
		r.toVector( result );

		ublas::subrange( jacobian, 0, 4, 0, 4 ) = ublas::identity_matrix< typename MT::value_type >( 4 );
	}

protected:
		double m_deltaTime;
		unsigned m_order;
};

} } } // namespace Ubitrack::Tracking::Function

#endif

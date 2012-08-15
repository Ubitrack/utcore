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
 * Function that performs a time update on a pose and its derivatives 
 * for the case of inside-out tracking of an object that is static
 * wrt. the world.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 


#ifndef __UBITRACK_TRACKING_FUNCTION_INSIDEOUTPOSETIMEUPDATE_H_INCLUDED__
#define __UBITRACK_TRACKING_FUNCTION_INSIDEOUTPOSETIMEUPDATE_H_INCLUDED__

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "LinearTimeUpdate.h"
#include "QuaternionTimeUpdate.h"
#include <utMath/Function/QuaternionVectorRotation.h>
#include <utMath/Function/RotationVelocityIntegration.h>

namespace Ubitrack { namespace Tracking { namespace Function {

/**
 * Updates a pose assuming constant velocity in both position and orientation.
 * We assume that a static object is observed by a moving (esp. rotating!) camera.
 * 
 * The state vector (p, dp, r, dr) is updated to (p', dp', r', dr') with
 *  - p'  = r quat(dr, dt) r*   p  (r quat(dr, dt) r*)* + dp dt
 *  - dp' = r quat(dr, dt) r*  dp  (r quat(dr, dt) r*)*  (optional)
 *  - r'  = r quat(dr, dt)
 *  - dr' = dr
 *  - dt  = time interval
 *  - quat( dr, dt ) = 3-vector to quaternion conversion
 */
class InsideOutPoseTimeUpdate
{
public:
	/**
	 * constructor.
	 * @param deltaTime time to forward in seconds
	 * @param posOrder can only be 0 or 1
	 */
	InsideOutPoseTimeUpdate( double deltaTime, int posOrder )
		: m_deltaTime( deltaTime )
		, m_posOrder( posOrder )
		, m_rotOrder( 1 )
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
		typedef typename VT1::value_type T;
		namespace ublas = boost::numeric::ublas;
		const int rs = 3 + 3 * m_posOrder; // start of rotation
		const int ts = size(); // total size
		
		// update position

		// rotate angular velocity into the translation coordinate frame
		ublas::vector< T > angVelRotated( 3 );
		ublas::matrix< T, ublas::column_major > jAngVelRotatedQ( 3, 4 );
		ublas::matrix< T, ublas::column_major > jAngVelRotatedV( 3, 3 );
		Math::Function::QuaternionVectorRotation().evaluateWithJacobian( angVelRotated,
			ublas::subrange( input, rs, rs + 4 ), ublas::subrange( input, rs + 4, rs + 7 ),
			jAngVelRotatedQ, jAngVelRotatedV );

		// integrate angular velocity
		ublas::vector< T > transUpdateRotation( 4 );
		ublas::matrix< T, ublas::column_major > jTransRotIntegrate( 4, 3 );
		Math::Function::RotationVelocityIntegration( m_deltaTime ).evaluateWithJacobian( 
			transUpdateRotation, angVelRotated, jTransRotIntegrate );

		// rotate translation by new quaternion
		ublas::vector< T > newTranslationTmp( 3 );
		ublas::matrix< T, ublas::column_major > jRotateTranslationQ( 3, 4 );
		ublas::matrix< T, ublas::column_major > jRotateTranslationV( 3, 3 );
		Math::Function::QuaternionVectorRotation().evaluateWithJacobian( newTranslationTmp,
			transUpdateRotation, ublas::subrange( input, 0, 3 ), jRotateTranslationQ, jRotateTranslationV );

		// add translation caused by constant velocity
		if ( m_posOrder == 1 )
			newTranslationTmp += m_deltaTime * ublas::subrange( input, 3, 6 );

		ublas::subrange( result, 0, 3 ) = newTranslationTmp;

		// jacobian d translation / d translation
		noalias( ublas::subrange( jacobian, 0, 3, 0, 3 ) ) = jRotateTranslationV;

		// jacobian d translation / d speed
		if ( m_posOrder == 1 )
			noalias( ublas::subrange( jacobian, 0, 3, 3, 6 ) ) = ublas::identity_matrix< T >( 3 ) * m_deltaTime;

		// jacobian d translation/d rotation
		{
		ublas::matrix< T, ublas::column_major > jAcc( ublas::prod( jRotateTranslationQ, jTransRotIntegrate ) );
		noalias( ublas::subrange( jacobian, 0, 3, rs, rs + 4 ) ) = ublas::prod( jAcc, jAngVelRotatedQ );

		// jacobian d translation / d angular velocity
		noalias( ublas::subrange( jacobian, 0, 3, rs + 4, rs + 7 ) ) = ublas::prod( jAcc, jAngVelRotatedV );
		}


		if ( m_posOrder == 1 )
		{
			// also rotate velocity
			Math::Function::QuaternionVectorRotation().evaluateWithJacobian( newTranslationTmp,
				transUpdateRotation, ublas::subrange( input, 3, 6 ), jRotateTranslationQ, jRotateTranslationV );
			ublas::subrange( result, 3, 6 ) = newTranslationTmp;

			// jacobian d velocity / d translation
			noalias( ublas::subrange( jacobian, 3, 6, 0, 3 ) ) = ublas::zero_matrix< T >( 3, 3 );

			// jacobian d velocity / d velocity
			noalias( ublas::subrange( jacobian, 3, 6, 3, 6 ) ) = jRotateTranslationV;

			// jacobian d velocity/d rotation
			ublas::matrix< T, ublas::column_major > jAcc( ublas::prod( jRotateTranslationQ, jTransRotIntegrate ) );
			noalias( ublas::subrange( jacobian, 3, 6, rs, rs + 4 ) ) = ublas::prod( jAcc, jAngVelRotatedQ );

			// jacobian d velocity / d angular velocity
			noalias( ublas::subrange( jacobian, 3, 6, rs + 4, rs + 7 ) ) = ublas::prod( jAcc, jAngVelRotatedV );
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

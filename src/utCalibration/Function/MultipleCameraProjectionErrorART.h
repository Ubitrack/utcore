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
 * @ingroup calibration
 * @file
 * functions for 3D->2D projections
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include "QuaternionRotationError.h"
#include "Dehomogenization.h"
 
#ifndef __UBITRACK_CALIBRATION_FUNCTION_MULTIPLECAMERAPROJECTIONERRORART_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_MULTIPLECAMERAPROJECTIONERRORART_H_INCLUDED__
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Jacobian for computing the pose error resulting from a projection of multiple 3d points with a single camera. 
 * Pose Error is expressed in A.R.T. format (which is different from Ubitrack)
 *
 * For each 3D point p, the jacobian of the projection
@verbatim
dehomogenize( P * ( E_r * R * ( p - c_g ) + R * c_g + t + e_t ) )
@endverbatim
 * is computed w.r.t (e_tx, e_ty, e_tz, e_rx, e_ry, e_rz) where (e_rx, e_ry, e_rz) are small rotation angles in around 
 * the x, y and z axes, given in radians.
 * E_r is assumed have expectation I and e_t 0. 
 * \c P is the 3x4 projection matrix of each camera, \c R the orientation (as an exponential map 3-vector), \c t the translation.
 *
 * p and P must be already known, the 6-vector (t, r) is the input to the function.
 */
template< class VType >
class MultipleCameraProjectionErrorART
{
public:
	/** 
	 * constructor.
	 * note: all parameters must stay constant during lifetime of the object
	 * @param p reference to vector of 3D-points to be projected 
	 * @param cameras reference to a vector of 3x4 camera matrices
	 * @param visibilities reference to a vector of observations. Each vector element contains a 
	 *    pair (i_p, i_c) which specifies that camera i_c has measured point i_p.
	 * @param centerOfGravity point in body coordinates that is used as the origin of the error
	 */
	MultipleCameraProjectionErrorART( const std::vector< Math::Vector< VType, 3 > >& p3D, 
		const std::vector< Math::Matrix< VType, 3, 4 > >& cameras, 
		const std::vector< std::pair< unsigned, unsigned > > visibilities,
		const Math::Vector< VType, 3 >& centerOfGravity = Math::Vector< VType, 3 >::zeros() )
		: m_p3D( p3D )
		, m_cam( cameras )
		, m_vis( visibilities )
		, m_centerOfGravity( centerOfGravity )
	{}

	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 2 * m_vis.size(); }

	/**
	 * @param param containing the pose as 6-vector, where rotation is stored in elements 4-6 as exponential map
	 * @param j matrix to store the jacobian (evaluated for param) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		using namespace Ubitrack::Math;
		namespace ublas = boost::numeric::ublas;

		// convert quaternion to matrix (for speedup)
		Quaternion rotQ( Quaternion::fromLogarithm( ublas::subrange( input, 3, 6 ) ) );
		Matrix< VType, 3, 3 > rot( rotQ );
		
		// create matrices
		Matrix< VType, 2, 3 > projJ;
		Matrix< VType, 3, 3 > rotJ;
		Vector< VType, 3 > rotated;
		Vector< VType, 3 > translated;
		Vector< VType, 3 > projected;
		
		for ( unsigned i = 0; i < m_vis.size(); i++ )
		{
			// shortcuts
			const Math::Vector< VType, 3 >& p3D( m_p3D[ m_vis[ i ].first ] );
			const Math::Matrix< VType, 3, 4 >& cam( m_cam[ m_vis[ i ].second ] );
			
			// rotate & project points
			noalias( rotated ) = ublas::prod( rot, p3D - m_centerOfGravity );
			noalias( translated ) = rotated + ublas::subrange( input, 0, 3 ) + ublas::prod( rot, m_centerOfGravity );
			noalias( projected ) = ublas::prod( ublas::subrange( cam, 0, 3, 0, 3 ), translated ) + ublas::column( cam, 3 );
			
			// create rotation jacobian for this measurement
			rotJ( 0, 0 ) = 0;
			rotJ( 0, 1 ) = rotated( 2 );
			rotJ( 0, 2 ) = -rotated( 1 );
			rotJ( 1, 0 ) = -rotated( 2 );
			rotJ( 1, 1 ) = 0;
			rotJ( 1, 2 ) = rotated( 0 );
			rotJ( 2, 0 ) = rotated( 1 );
			rotJ( 2, 1 ) = -rotated( 0 );
			rotJ( 2, 2 ) = 0;
			
			// create projection jacobian for this measurement
			Dehomogenization< 3 >().jacobian( projected, projJ );
			
			// apply chain rule
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 0, 3 ) ) 
				= ublas::prod( projJ, ublas::subrange( cam, 0, 3, 0, 3 ) );
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 3, 6 ) ) 
				= ublas::prod( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 0, 3 ), rotJ );
		}
	}
	
protected:
	const std::vector< Math::Vector< VType, 3 > >& m_p3D;
	const std::vector< Math::Matrix< VType, 3, 4 > >& m_cam;
	const std::vector< std::pair< unsigned, unsigned > > m_vis;
	const Math::Vector< VType, 3 > m_centerOfGravity;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

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
 
#ifndef __UBITRACK_CALIBRATION_FUNCTION_MULTIPLECAMERAPROJECTIONERROR_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_MULTIPLECAMERAPROJECTIONERROR_H_INCLUDED__
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Jacobian for computing the pose error resulting from a projection of multiple 3d points with a single camera. 
 *
 * For each 3D point p, the jacobian of the projection
@verbatim
dehomogenize( C * ( r * e_r * p * e_r' * r' + t + e_t ) )
@endverbatim
 * is computed w.r.t (e_tx, e_ty, e_tz, e_rx, e_ry, e_rz) where (e_rx, e_ry, e_rz) is the imaginary part of e_r.
 * e_r is assumed to be (0, 0, 0, 1) and e_t is 0. 
 * \c c is the 3x4 projection matrix of each camera, \c r the orientation (as a quaternion), \c t the translation.
 *
 * p and C must be already known, the 7-vector (t, r) is the input to the function.
 */
template< class VType >
class MultipleCameraProjectionError
{
public:
	/** 
	 * constructor.
	 * note: all parameters must stay constant during lifetime of the object
	 * @param p reference to vector of 3D-points to be projected 
	 * @param cameras reference to a vector of 3x4 camera matrices
	 * @param visibilities reference to a vector of observations. Each vector element contains a 
	 *    pair (i_p, i_c) which specifies that camera i_c has measured point i_p.
	 */
	MultipleCameraProjectionError( const std::vector< Math::Vector< 3, VType > >& p3D, 
		const std::vector< Math::Matrix< 3, 4, VType > >& cameras, 
		const std::vector< std::pair< std::size_t, std::size_t > > visibilities )
		: m_p3D( p3D )
		, m_cam( cameras )
		, m_vis( visibilities )
	{}

	/**
	 * return the size of the result vector
	 */
	std::size_t size() const
	{ return 2 * m_vis.size(); }

	/**
	 * @param param containing the parameters (pose as 7-vector>
	 * @param j matrix to store the jacobian (evaluated for param) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		using namespace Ubitrack::Math;
		namespace ublas = boost::numeric::ublas;

		// convert quaternion to matrix (for speedup)
		Quaternion rotQ( Quaternion::fromVector( ublas::subrange( input, 3, 7 ) ) );
		Matrix< 3, 3, VType > rot( rotQ );
		
		// create matrices
		Matrix< 2, 3, VType > projJ;
		Matrix< 3, 3, VType > rotJ;
		Vector< 3, VType > rotated;
		Vector< 3, VType > projected;
		
		for ( std::size_t i( 0 ); i < m_vis.size(); i++ )
		{
			// shortcuts
			const Math::Vector< 3, VType >& p3D( m_p3D[ m_vis[ i ].first ] );
			const Math::Matrix< 3, 4, VType >& cam( m_cam[ m_vis[ i ].second ] );
			
			// rotate & project points
			noalias( rotated ) = ublas::prod( rot, p3D ) + ublas::subrange( input, 0, 3 );
			noalias( projected ) = ublas::prod( ublas::subrange( cam, 0, 3, 0, 3 ), rotated ) + ublas::column( cam, 3 );
			
			// create jacobian for this measurement
			QuaternionRotationError< VType > qrf( p3D );
			qrf.jacobian( ublas::subrange( input, 3, 7 ), rotJ );
			Dehomogenization< 3 >().jacobian( projected, projJ );
			
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 0, 3 ) ) 
				= ublas::prod( projJ, ublas::subrange( cam, 0, 3, 0, 3 ) );
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 3, 6 ) ) 
				= ublas::prod( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 0, 3 ), rotJ );
		}
	}
	
protected:
	const std::vector< Math::Vector< 3, VType > >& m_p3D;
	const std::vector< Math::Matrix< 3, 4, VType > >& m_cam;
	const std::vector< std::pair< std::size_t, std::size_t > > m_vis;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

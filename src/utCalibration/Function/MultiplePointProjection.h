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

#ifndef __UBITRACK_CALIBRATION_FUNCTION_MULTIPLEPOINTPROJECTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_MULTIPLEPOINTPROJECTION_H_INCLUDED__
 
#include "QuaternionRotation.h"
#include "Dehomogenization.h"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * Function that projects multiple 3D points to 2D points given the pose. 
 *
 * It computes for each point in the list 
@verbatim
dehomogenize( C * ( r * p * r' + t ) )
@endverbatim
 * and/or the jacobian of this function with respect to (t, r), where \c C is the camera 
 * intrinsics matrix, \c r the orientation (as a quaternion), \c t the translation and
 * \c p the 3D point to be projected.
 *
 * p and C must be already known, the 7-vector (t, r) is the input to the function.
 *
 * This function is used in Pose Estimation and error propagation.
 */
template< class VType >
class MultiplePointProjection
{
public:
	/** 
	 * constructor.
	 * @param p reference to vector of 3D-points to be projected (must stay constant during lifetime of the object)
	 * @param cam reference to 3x3 camera intrinsics matrix (must stay constant during lifetime of the object)
	 */
	MultiplePointProjection( const std::vector< Math::Vector< 3, VType > >& p3D, const Math::Matrix< 3, 3, VType >& cam )
		: m_p3D( p3D )
		, m_cam( cam )
	{}

	/**
	 * return the size of the result vector
	 */
	std::size_t size() const
	{ return 2 * m_p3D.size(); }

	/**
	 * @param result 2*N-vector to store the result in
	 * @param input containing the parameters (tx, ty, tz, qx, qy, qz, qw)
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		using namespace Ubitrack::Math;
		namespace ublas = boost::numeric::ublas;

		// convert quaternion to matrix (for speedup)
		Quaternion rotQ( Quaternion::fromVector( ublas::subrange( input, 3, 7 ) ) );
		Matrix< 3, 3, VType > rot( rotQ );
		
		for ( std::size_t i ( 0 ); i < m_p3D.size(); ++i )
		{
			// rotate & project points
			Vector< 3, VType > rotated( ublas::prod( rot, m_p3D[ i ] ) + ublas::subrange( input, 0, 3 ) );
			Vector< 3, VType > projected( ublas::prod( m_cam, rotated ) );
			ublas::noalias( ublas::subrange( result, i * 2, (i+1) * 2 ) ) = ublas::subrange( projected, 0, 2 ) / projected( 2 );
		}
	}
	
	/**
	 * @param result vector to store the result in
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		// TODO: implement as one function (more efficient)
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the jacobian (evaluated for input) in
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
		Matrix< 3, 4, VType > rotJ;
		Vector< 3, VType > rotated;
		Vector< 3, VType > projected;
		
		for ( unsigned i = 0; i < m_p3D.size(); i++ )
		{
			// rotate & project points
			noalias( rotated ) = ublas::prod( rot, m_p3D[ i ] ) + ublas::subrange( input, 0, 3 );
			noalias( projected ) = ublas::prod( m_cam, rotated );
			
			// create jacobian for this measurement
			QuaternionRotation< VType > qrf( m_p3D[ i ] );
			qrf.jacobian( ublas::subrange( input, 3, 7 ), rotJ );
			Dehomogenization< 3 >().jacobian( projected, projJ );
			
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 0, 3 ) ) 
				= ublas::prod( projJ, m_cam );
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 3, 7 ) ) 
				= ublas::prod( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 0, 3 ), rotJ );
		}
	}
	
protected:
	const std::vector< Math::Vector< 3, VType > >& m_p3D;
	const Math::Matrix< 3, 3, VType >& m_cam;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

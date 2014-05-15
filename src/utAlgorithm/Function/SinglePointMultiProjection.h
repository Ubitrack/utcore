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
 * functions for (n*2D)->3D estimations (n>1).
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_SINGLEPOINTPROJECTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_SINGLEPOINTPROJECTION_H_INCLUDED__
 
#include <utMath/Geometry/PointProjection.h>
#include <utMath/Geometry/PointTransformation.h>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace Ubitrack { namespace Algorithm { namespace Function {

/**
 * Function that projects a single 3D point using several projections.
 *
 * It computes for each 3x4 - projection in the container 
@verbatim
projection( [ R_i, T_i] * [p_x, p_y, p_z, 1]^t )
@endverbatim
 * and/or the Jacobian of this function with respect to [p_x, p_y, p_z] , where \c [R_i, T_i] is 
 * an extrinsic 3x4camera matrix, \c R the orientation , \c T the translation.
 *
 * R and T must be already known, the 3-vector [p_x, p_y, p_z] is the input to the function.
 *
 * This function is used in 3DPointReconstruction.
 */
template< class VType, typename ForwardIterator1 >
class SinglePointMultiProjection
{
protected:
	const ForwardIterator1 m_iBegin;
	const ForwardIterator1 m_iEnd;
	
public:
	/** 
	 * constructor.
	 * @param iBegin iterator to the beginning of a contianer with projections(must stay constant during lifetime of the object)
	 * @param iEnd iterator to the end of a container with projections(must stay constant during lifetime of the object)
	 */
	SinglePointMultiProjection( const ForwardIterator1 iBegin, const ForwardIterator1 iEnd )
		: m_iBegin( iBegin )
		, m_iEnd( iEnd )
	{}

	/**
	 * return the size of the result vector containing the reprojections
	 */
	std::size_t size() const
	{ return 2 * ( std::distance( m_iBegin, m_iEnd ) ); }

	/**
	 * @param result 2*N-vector to store the result in
	 * @param input containing the parameters ( px, py, pz )
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		namespace ublas = boost::numeric::ublas;

		std::size_t i ( 0 );
		for ( ForwardIterator1 it ( m_iBegin ); it != m_iEnd; ++i, ++it )
			ublas::subrange ( result, i*2+0, i*2+2 ) = Math::Geometry::ProjectPoint()( *it, static_cast< Math::Vector< VType, 3 > > ( ublas::subrange( input, 0, 3 )  ) );
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
		const Math::Vector< VType, 3 > vec( input( 0 ), input( 1 ), input( 2 ) );		
		std::size_t i( 0 );
		for ( ForwardIterator1 it ( m_iBegin ); it != m_iEnd; ++i, ++it )
		{
			
			const Math::Vector< VType, 3 > point ( Math::Geometry::TransformPoint()( *it, vec ) );
			const VType P3_14p = point( 2 );
			//= 1/(P3_14p^2)
			const VType P3_14p2 = 1 / (P3_14p * P3_14p); 
			const VType P1_14p = P3_14p2 * point( 0 );
			const VType P2_14p = P3_14p2 * point( 1 );
			
			const VType P3_1 = (*it)( 2, 0 );
			const VType P3_2 = (*it)( 2, 1 );
			const VType P3_3 = (*it)( 2, 2 );
			J( i*2+0, 0 ) = (*it)( 0, 0 ) / P3_14p - P3_1 * P1_14p;
			J( i*2+0, 1 ) = (*it)( 0, 1 ) / P3_14p - P3_2 * P1_14p;
			J( i*2+0, 2 ) = (*it)( 0, 2 ) / P3_14p - P3_3 * P1_14p;
			J( i*2+1, 0 ) = (*it)( 1, 0 ) / P3_14p - P3_1 * P2_14p;
			J( i*2+1, 1 ) = (*it)( 1, 1 ) / P3_14p - P3_2 * P2_14p;
			J( i*2+1, 2 ) = (*it)( 1, 2 ) / P3_14p - P3_3 * P2_14p;
		}
	}
};

} } } // namespace Ubitrack::Algorithm::Function

#endif //__UBITRACK_CALIBRATION_FUNCTION_SINGLEPOINTPROJECTION_H_INCLUDED__

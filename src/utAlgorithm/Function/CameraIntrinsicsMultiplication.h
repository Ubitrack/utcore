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
 * Multiplies a 2-vector with a camera intrinsics matrix given as a 5-vector
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_CAMERAINTRINSICSMULTIPLICATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_CAMERAINTRINSICSMULTIPLICATION_H_INCLUDED__
 
namespace Ubitrack { namespace Algorithm { namespace Function {

/**
 * Applies an intrinsic camera matrix (given as a 5-vector) to an already dehomogenized 2-vector.
 * The jacobian is computed wrt. the 5-vector representation of the matrix.
 * We assume that the lower right matrix element is -1. See src/Vision/Readme.txt why.
 */
template< typename T >
class CameraIntrinsicsMultiplication
{
public:
	/**
	 * Constructor.
	 * @param p reference to the point being transformed
	 */
	CameraIntrinsicsMultiplication( const Math::Vector< T, 2 >& p )
		: m_p( p )
	{}

	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 2; }

	/*
	 * @param result the transformed point (2-vector)
	 * @param input the 5-vector camera parameters 
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		result( 0 ) = - ( input( 0 ) * m_p( 0 ) + input( 1 ) * m_p( 1 ) + input( 2 ) );
		result( 1 ) = - (                         input( 3 ) * m_p( 1 ) + input( 4 ) );
	}
	
	/*
	 * @param result the transformed point (2-vector)
	 * @param input the 5-vector camera parameters 
	 * @param J jacobian (2x5-matrix)
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		evaluate( result, input );
		jacobian( input, J );
	}

	/*
	 * @param input the 5-vector camera parameters 
	 * @param J jacobian (2x5-matrix)
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		J( 0, 0 ) = -m_p( 0 );
		J( 0, 1 ) = -m_p( 1 );
		J( 0, 2 ) = -1;
		J( 0, 3 ) = 0;
		J( 0, 4 ) = 0;
		J( 1, 0 ) = 0;
		J( 1, 1 ) = 0;
		J( 1, 2 ) = 0;
		J( 1, 3 ) = -m_p( 1 );
		J( 1, 4 ) = -1;
	}
	
protected:
	const Math::Vector< T, 2 >& m_p;
};

} } } // namespace Ubitrack::Algorithm::Function

#endif

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
 * Radially distorts a 2-vector
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_RADIALDISTORTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_RADIALDISTORTION_H_INCLUDED__
 
namespace Ubitrack { namespace Calibration { namespace Function {

/**
 * \internal
 */
template< class VT1, class VT2, class VT3 >
void radialDistortion( VT1& result, const VT2& p, const VT3& d )
{
	typename VT1::value_type r2 = p( 0 ) * p( 0 ) + p( 1 ) * p( 1 );
	noalias( result ) = p * ( 1 + d( 0 ) * r2 + d( 1 ) * r2 * r2 );
	result( 0 ) +=  2 * d ( 2 ) * p( 0 ) * p( 1 ) + d( 3 ) * ( r2 + 2 * p( 0 ) * p( 0 ) );
	result( 1 ) +=  2 * d ( 3 ) * p( 0 ) * p( 1 ) + d( 2 ) * ( r2 + 2 * p( 1 ) * p( 1 ) );
}


/**
 * Radially and tangentially distorts a 2-vector (x,y). The distortion is described by a 4-vector (k1, k2, p1, p2):
 *    x' = x + x( k1 * r^2 + k2 * r^4 ) + ( 2 * p1 * x * y + p2 * ( r^2 + 2 * x^2 ) )
 *    y' = y + y( k1 * r^2 + k2 * r^4 ) + ( 2 * p2 * x * y + p1 * ( r^2 + 2 * y^2 ) )
 * where
 *    r^2 = x^2 + y^2
 *
 * The jacobian is computed wrt. the point p.
 */
template< typename T >
class RadialDistortionWrtP
{
public:
	/**
	 * Constructor.
	 * @param d reference to the distortion coefficients
	 */
	RadialDistortionWrtP( const Math::Vector< T, 4 >& d )
		: m_d( d )
	{}

	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 2; }

	/**
	 * @param result the transformed point (2-vector)
	 * @param input the 2-vector describing the undistorted point
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{ radialDistortion( result, input, m_d ); }
	
	/**
	 * @param result the transformed point (2-vector)
	 * @param input the 2-vector describing the undistorted point
	 * @param J jacobian (2x2-matrix)
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input the 2-vector describing the undistorted point
	 * @param J jacobian (2x2-matrix)
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		typedef typename MT::value_type VType;
		
		VType t331 = input( 0 ) * input( 0 );
		VType t332 = input( 1 ) * input( 1 );
		VType t340 = m_d( 2 ) * input( 0 );
		VType t341 = t331 + t332;
		VType t342 = 2 * m_d( 1 ) * t341;
		VType t343 = m_d( 0 ) + t342;
		VType t344 = input( 0 ) * t343;
		VType t345 = m_d( 3 ) + t344;
		VType t346 = input( 1 ) * t345;
		VType t347 = t340 + t346;
		VType t348 = 2 * t347;
		VType t328 = t331 * t331;
		VType t333 = 6 * m_d( 1 ) * t331 * t332;
		VType t334 = t332 * t332;
		J( 0, 0 ) = 1 + 6 * m_d( 3 ) * input( 0 ) + 2 * m_d( 2 ) * input( 1 ) + 5 * m_d( 1 ) * t328 +
			m_d( 0 ) * ( 3 * t331 + t332 ) + t333 + m_d( 1 ) * t334;
		J( 0, 1 ) = t348;
		J( 1, 0 ) = t348;
		J( 1, 1 ) = 1 + 2 * m_d( 3 ) * input( 0 ) + 6 * m_d( 2 ) * input( 1 ) + m_d( 1 ) * t328 + 
			m_d( 0 ) * ( t331 + 3 * t332 ) + t333 + 5 * m_d( 1 ) * t334;	
	}
	
protected:
	const Math::Vector< T, 4 >& m_d;
};


/**
 * Radially and tangentially distorts a 2-vector (x,y). The distortion is described by a 4-vector (k1, k2, p1, p2):
 *    x' = x + x( k1 * r^2 + k2 * r^4 ) + ( 2 * p1 * x * y + p2 * ( r^2 + 2 * x^2 ) )
 *    y' = y + y( k1 * r^2 + k2 * r^4 ) + ( 2 * p2 * x * y + p1 * ( r^2 + 2 * y^2 ) )
 * where
 *    r^2 = x^2 + y^2
 *
 * The jacobian is computed wrt. the distortion parameters d.
 */
template< typename T >
class RadialDistortionWrtD
{
public:
	/**
	 * Constructor.
	 * @param p reference to the distortion coefficients
	 */
	RadialDistortionWrtD( const Math::Vector< T, 2 >& p )
		: m_p( p )
	{}

	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 2; }

	/**
	 * @param result the transformed point (2-vector)
	 * @param input the 4-vector of distortion parameters
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{ radialDistortion( result, m_p, input ); }
	
	/**
	 * @param result the transformed point (2-vector)
	 * @param input the 4-vector of distortion parameters
	 * @param J jacobian (2x4-matrix)
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input the 4-vector of distortion parameters
	 * @param J jacobian (2x2-matrix)
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		typedef typename MT::value_type VType;
		
		VType t380 = m_p( 0 ) * m_p( 0 );
		VType t381 = m_p( 1 ) * m_p( 1 );
		VType t382 = t380 + t381;
		VType t384 = t382 * t382;
		VType t386 = 2 * m_p( 0 ) * m_p( 1 );
		J( 0, 0 ) = m_p( 0 ) * t382;
		J( 0, 1 ) = m_p( 0 ) * t384;
		J( 0, 2 ) = t386;
		J( 0, 3 ) = 3 * t380 + t381;
		J( 1, 0 ) = m_p( 1 ) * t382;
		J( 1, 1 ) = m_p( 1 ) * t384;
		J( 1, 2 ) = t380 + 3 * t381;
		J( 1, 3 ) = t386;
	}
	
protected:
	const Math::Vector< T, 2 >& m_p;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

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
 * defines a linear function used in rotation hand-eye-calibration
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
namespace Ubitrack { namespace Algorithm { namespace Function {

/**
 * Defines the function used in optimization of rotation-only hand-eye-calibration.
 * Computes a x - x b (which should be 0) and its jacobian.
 * Modeled after UnaryFunctionPrototype.
 */
class RotHecMeasurement
{
public:
	/** construct the object, given the a and b measurements (which are not being optimized) */
	RotHecMeasurement( const Math::Quaternion& a, const Math::Quaternion& b )
		: m_a( a ), m_b( b )
	{}
	
	unsigned size() const
	{ return 4; }

	/**
	 * @param result 4-vector, stores the result of ax-xb
	 * @param input 4-vector, containing x as a quaternion (x, y, z, w)
	 * @param jacobian 4x4-matrix where the resulting jacobian is stored
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& jacobian ) const
	{
		// compute result
		Math::Quaternion x( Math::Quaternion::fromVector( input ) );
		Math::Quaternion z = m_a * x - x * m_b;
		z.toVector( result );
		
		// compute jacobian
		jacobian( 0, 0 ) =  m_a.w() - m_b.w();
		jacobian( 0, 1 ) = -m_a.z() - m_b.z();
		jacobian( 0, 2 ) =  m_a.y() + m_b.y();
		jacobian( 0, 3 ) =  m_a.x() - m_b.x();

		jacobian( 1, 0 ) =  m_a.z() + m_b.z();
		jacobian( 1, 1 ) =  m_a.w() - m_b.w();
		jacobian( 1, 2 ) = -m_a.x() - m_b.x();
		jacobian( 1, 3 ) =  m_a.y() - m_b.y();

		jacobian( 2, 0 ) = -m_a.y() - m_b.y();
		jacobian( 2, 1 ) =  m_a.x() + m_b.x();
		jacobian( 2, 2 ) =  m_a.w() - m_b.w();
		jacobian( 2, 3 ) =  m_a.z() - m_b.z();

		jacobian( 3, 0 ) = -m_a.x() + m_b.x();
		jacobian( 3, 1 ) = -m_a.y() + m_b.y();
		jacobian( 3, 2 ) = -m_a.z() + m_b.z();
		jacobian( 3, 3 ) =  m_a.w() - m_b.w();
	}
	
protected:
	Math::Quaternion m_a;
	Math::Quaternion m_b;
};

/**
 * used to compute the covariance of a x - x b, given x and covariances of a and b
 */
class RotHecCombine
{
public:
	RotHecCombine( const Math::Quaternion& x )
		: m_x( x )
	{}
	
	unsigned size() const
	{ return 4; }

	/**
	 * @param jacobian1 the jacobian of ax-xb wrt. a (output, 4x4-matrix)
	 * @param jacobian2 the jacobian of ax-xb wrt. b (output, 4x4-matrix)
	 */
	template< class VT2, class VT3, class MT1, class MT2 > 
	void jacobian( const VT2&, const VT3&, MT1& jacobian1, MT2& jacobian2 ) const
	{
		jacobian1( 0, 0 ) =  m_x.w();
		jacobian1( 0, 1 ) =  m_x.z();
		jacobian1( 0, 2 ) = -m_x.y();
		jacobian1( 0, 3 ) =  m_x.x();
		
		jacobian1( 1, 0 ) = -m_x.z();
		jacobian1( 1, 1 ) =  m_x.w();
		jacobian1( 1, 2 ) =  m_x.x();
		jacobian1( 1, 3 ) =  m_x.y();

		jacobian1( 2, 0 ) =  m_x.y();
		jacobian1( 2, 1 ) = -m_x.x();
		jacobian1( 2, 2 ) =  m_x.w();
		jacobian1( 2, 3 ) =  m_x.z();

		jacobian1( 3, 0 ) = -m_x.x();
		jacobian1( 3, 1 ) = -m_x.y();
		jacobian1( 3, 2 ) = -m_x.z();
		jacobian1( 3, 3 ) =  m_x.w();
		
		jacobian2 = -jacobian1;
	}
	
protected:
	Math::Quaternion m_x;
};

} } } // namespace Ubitrack::Algorithm::Function

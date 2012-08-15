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
 * @ingroup math
 * @file
 * function class for integrating a rotation velocity 3-vector
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_ROTATIONVELOCITYINTEGRATION_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_ROTATIONVELOCITYINTEGRATION_H_INCLUDED__

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace Ubitrack { namespace Math { namespace Function {
 
/**
 * Function class for integrating a rotation velocity 3-vector.
 * Note: The computation of the jacobian assumes a "small" result!
 *
 */
struct RotationVelocityIntegration
{
	/**
	 * Constructor.
	 * @param dt the time to integrate
	 */
	RotationVelocityIntegration( double dt )
		: m_dt( dt )
	{}

	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 4; }
	
	/**
	 * Evaluate the function on the input \c input and store the result in
	 * \c result. \c VT1, \c VT2 and \c VT3 can be assumed to be Boost uBlas vectors
	 * using operator() to access elements
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& q ) const
	{
		typedef typename VT1::value_type T;
		T fNorm = boost::numeric::ublas::norm_2( q );
		T s;
		if ( fabs( fNorm * m_dt ) > 1e12 )
			s = sin( fNorm * m_dt * T( 0.5 ) ) / fNorm;
		else
			s = m_dt * T( 0.5 );
		boost::numeric::ublas::subrange( result, 0, 3 ) = q * s;
		result( 3 ) = cos( fNorm * m_dt * T( 0.5 ) );
	}
	
	/**
	 * Evaluate the function on the input \c input and return both the result
	 * and the jacobian. All parameters can be assumed to be Boost 
	 * uBlas vectors and matrices using operator() to access elements.
	 */
	template< class VT1, class VT2, class MT1 > 
	void evaluateWithJacobian( VT1& result, const VT2& v, MT1& jac ) const
	{
		jacobian( v, jac );
		evaluate( result, v );
	}
	
	/**
	 * Evaluate the jacobian on the input \c input.
	 * All parameters can be assumed to be Boost uBlas vectors and matrices 
	 * using operator() to access elements.
	 */
	template< class VT2, class MT1 > 
	void jacobian( const VT2& v, MT1& jac ) const
	{
		typedef typename MT1::value_type T;
		boost::numeric::ublas::subrange( jac, 0, 3, 0, 3 ) = T( m_dt * 0.5 ) * boost::numeric::ublas::identity_matrix< T >( 3 );
		boost::numeric::ublas::row( jac, 3 ) = T( m_dt * -0.25 ) * v;
	}

	double m_dt;
};

} } } // namespace Ubitrack::Math::Function

#endif


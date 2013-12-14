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
 * Class that numerically approximates the jacobian of a function.
 *
 * This class will automatically add the evaluateWithJacobian() and jacobian() methods
 * to any function class that only implements evaluate().
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include <utMath/Vector.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
 
namespace Ubitrack { namespace Math { namespace Function {
 
/**
 * Function class that numerically approximates the jacobian of a function.
 * This requires n function evaluations for each jacobian computation, where n is the size of the
 * input vector.
 */
template< class FC >
class DiscreteJacobianApproximation
{
public:
	/**
	 * construct a new approximation.
	 * @param f the function object whose jacobian is to be estimated
	 * @param fApproxWith distance of the differencing points expressed as a percentage 
	 * of the absolute value of each parameter
	 */
	DiscreteJacobianApproximation( const FC& f, double fApproxWidth = 0.001 )
		: m_f( f )
		, m_fApproxWidth( fApproxWidth )
	{}
	
	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return m_f.size(); }
	

	/**
	 * Evaluate the function on the input \c input and store the result in
	 * \c result. \c VT1 and \c VT2 can be assumed to be Boost uBlas vectors
	 * using operator() to access elements
	 *
	 * @param result vector to store the result in
	 * @param input containing the parameters (to be optimized)
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{ return m_f.evaluate( result, input ); }
	
	/**
	 * Evaluate the function on the input \c input and return both the resultscons v
	 * and the jacobian. \c VT1, \c VT2 and \c MT can be assumed to be Boost 
	 * uBlas vectors and matrices using operator() to access elements.
	 *
	 * @param result vector to store the result in
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		namespace ublas = boost::numeric::ublas;
		
		// evaluate at initial position
		m_f.evaluate( result, input );
		
		// evaluate for each discretization point
		unsigned inputSize = input.size();
		ublas::vector< typename VT1::value_type > testResult( m_f.size() );
		ublas::vector< typename VT2::value_type > testInput( input );
		for ( unsigned i = 0; i < inputSize; i++ )
		{
			// slightly modify the input
			double eps;
			if ( testInput( i ) != 0 )
				eps = testInput( i ) * m_fApproxWidth;
			else
				eps = m_fApproxWidth;

			testInput( i ) = input( i ) + eps;
				
			m_f.evaluate( testResult, testInput );
			
			// compute jacobian column
			ublas::column( J, i ) = ( testResult - result ) / eps;
			
			// reset testInput
			testInput( i ) = input( i );
		}
	}

	/**
	 * Compute only the jacobian evaluated at the given state.
	 *
	 * This is usually used in error propagation.
	 *
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		boost::numeric::ublas::vector< typename VT2::value_type > result( m_f.size() );
		evaluateWithJacobian( result, input, J );
	}
	
protected:
	FC m_f;
	double m_fApproxWidth;
};

} } } // namespace Ubitrack::Math::Function

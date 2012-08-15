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
 * prototype class for unary functions with derivatives
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

namespace Ubitrack { namespace Math { namespace Function {
 
/**
 * Prototype for unary functions used for optimization and error propagation.
 *
 * This class is not meant to be used directly, but shows the interface every 
 * unary function implementation has to provide. 
 *
 * You do not need to explicitly derive from this class!
 *
 * In most cases, implementing only the \c evaluateWithJacobian is sufficient.
 */
struct UnaryFunctionPrototype
{
	/**
	 * return the size of the result vector
	 */
	unsigned size() const;

	/**
	 * Evaluate the function on the input \c input and store the result in
	 * \c result. \c VT1 and \c VT2 can be assumed to be Boost uBlas vectors
	 * using operator() to access elements
	 *
	 * @param result vector to store the result in
	 * @param input containing the parameters (to be optimized)
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const;
	
	/**
	 * Evaluate the function on the input \c input and return both the result
	 * and the jacobian. \c VT1, \c VT2 and \c MT can be assumed to be Boost 
	 * uBlas vectors and matrices using operator() to access elements.
	 *
	 * @param result vector to store the result in
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const;

	/**
	 * Compute only the jacobian evaluated at the given state.
	 *
	 * This is usually used in error propagation.
	 *
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const;
};

} } } // namespace Ubitrack::Math::Function

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
 * prototype class for binary functions with derivatives
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

namespace Ubitrack { namespace Math { namespace Function {
 
/**
 * Prototype for binary functions used for error propagation.
 *
 * This class is not meant to be used directly, but shows the interface every 
 * unary function implementation has to provide. 
 *
 * You do not need to explicitly derive from this class!
 *
 * In most cases, implementing only the \c evaluateWithJacobian is sufficient.
 */
struct BinaryFunctionPrototype
{
	/**
	 * return the size of the result vector
	 */
	unsigned size() const;
	
	/**
	 * Evaluate the function on the input \c input and store the result in
	 * \c result. \c VT1, \c VT2 and \c VT3 can be assumed to be Boost uBlas vectors
	 * using operator() to access elements
	 *
	 * @param result vector to store the result in
	 * @param input1 vector containing the parameters of the first input
	 * @param input2 vector containing the parameters of the second input
	 */
	template< class VT1, class VT2, class VT3 > 
	void evaluate( VT1& result, const VT2& input1, const VT3& input2 ) const;
	
	/**
	 * Evaluate the function on the input \c input and return both the result
	 * and the jacobian. All parameters can be assumed to be Boost 
	 * uBlas vectors and matrices using operator() to access elements.
	 *
	 * @param result vector to store the result in
	 * @param input1 contains the first parameter vector
	 * @param input2 contains the second parameter vector
	 * @param jacobian1 matrix to store the jacobian wrt. the first vector in
	 * @param jacobian2 matrix to store the jacobian wrt. the second vector in
	 */
	template< class VT1, class VT2, class VT3, class MT1, class MT2 > 
	void evaluateWithJacobian( VT1& result, const VT2& input1, const VT3& input2, MT1& jacobian1, MT2& jacobian2 ) const;
	
	/**
	 * Evaluate the jacobian on the input \c input.
	 * All parameters can be assumed to be Boost uBlas vectors and matrices 
	 * using operator() to access elements.
	 *
	 * @param input1 contains the first parameter vector
	 * @param input2 contains the second parameter vector
	 * @param jacobian1 matrix to store the jacobian wrt. the first vector in
	 * @param jacobian2 matrix to store the jacobian wrt. the second vector in
	 */
	template< class VT2, class VT3, class MT1, class MT2 > 
	void jacobian( const VT2& input1, const VT3& input2, MT1& jacobian1, MT2& jacobian2 ) const;
};

} } } // namespace Ubitrack::Math::Function

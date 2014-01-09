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
 * @ingroup Math
 * @file
 * Class to derive multivariate functions from
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_MULTIVARIATEFUNCTION_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_MULTIVARIATEFUNCTION_H_INCLUDED__

#include "Detail/Binder.h"

namespace Ubitrack { namespace Math { namespace Optimization { namespace Function {

/**
 * Derive your own multivariate functions from this class.
 * @param Derived the class of the derived function
 * @param Size size of the result vector. 0 if unknown at compile-time.
 */
template< class Derived, unsigned Size >
class MultiVariateFunction
{
public:
	unsigned size() const
	{ return staticSize; }

private:
	template< class, class > friend class Detail::Binder;
	
	static const unsigned staticSize = Size;
	static const bool wantsJacobian = false;
	
	// functions to strip off parameter vector in evaluations

	// one-parameter function
	template< class ParameterVector, class DestinationVector, class Param1 >
	void i_evaluate( const ParameterVector& p, DestinationVector& d, const Param1& p1 ) const
	{
		// just strip parameter vector
		static_cast< const Derived* >( this )->evaluate( d, p1 );
	}
	
	// two-parameter function
	template< class ParameterVector, class DestinationVector, class Param1, class Param2 >
	void i_evaluate( const ParameterVector& p, DestinationVector& d, const Param1& p1, const Param2& p2 ) const
	{
		// just strip parameter vector
		static_cast< const Derived* >( this )->evaluate( d, p1, p2 );
	}

	// three-parameter function
	template< class ParameterVector, class DestinationVector, class Param1, class Param2, class Param3 >
	void i_evaluate( const ParameterVector&, DestinationVector& d, const Param1& p1, const Param2& p2, const Param3& p3 ) const
	{
		// just strip parameter vector
		static_cast< const Derived* >( this )->evaluate( d, p1, p2, p3 );
	}

	// stop recursion for: only evaluate parameters
	template< class ParameterVector >
	void i_evaluateParameters( const ParameterVector&  ) const
	{}

	
	// functions to strip off parameter vector in jacobian calculations

	// one-parameter function
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1 >
	void i_multiplyJacobian1( const ParameterVector&, const LeftHand& l, DestinationMatrix& j, const Param1& p1 ) const
	{
		static_cast< const Derived* >( this )->multiplyJacobian1( l, j, p1 );
	}

	// two-parameter functions
	
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void i_multiplyJacobian1( const ParameterVector&, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2 ) const
	{
		static_cast< const Derived* >( this )->multiplyJacobian1( l, j, p1, p2 );
	}
	
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void i_multiplyJacobian2( const ParameterVector&, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2 ) const
	{
		static_cast< const Derived* >( this )->multiplyJacobian2( l, j, p1, p2 );
	}

	// three-parameter functions
	
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2, class Param3 >
	void i_multiplyJacobian1( const ParameterVector&, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2, const Param3& p3 ) const
	{
		static_cast< const Derived* >( this )->multiplyJacobian1( l, j, p1, p2, p3 );
	}
	
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2, class Param3 >
	void i_multiplyJacobian2( const ParameterVector&, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2, const Param3& p3 ) const
	{
		static_cast< const Derived* >( this )->multiplyJacobian2( l, j, p1, p2, p3 );
	}
	
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2, class Param3 >
	void i_multiplyJacobian3( const ParameterVector&, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2, const Param3& p3 ) const
	{
		static_cast< const Derived* >( this )->multiplyJacobian3( l, j, p1, p2, p3 );
	}

	
	// terminate chain of calls to i_multiplyJacobian
	template< unsigned LHSize, class ParameterVector, class LeftHand, class DestinationMatrix, class Param1 >
	void i_multiplyJacobian( const ParameterVector&, const LeftHand&, DestinationMatrix&, const Param1&  ) const
	{}
	
	template< unsigned LHSize, class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void i_multiplyJacobian( const ParameterVector&, const LeftHand&, DestinationMatrix&, const Param1&, const Param2&  ) const
	{}
	
	template< unsigned LHSize, class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2, class Param3 >
	void i_multiplyJacobian( const ParameterVector&, const LeftHand&, DestinationMatrix&, const Param1&, const Param2&, const Param3&  ) const
	{}
};


}}}} // namespace Ubitrack::Math::Optimization::Function

#endif

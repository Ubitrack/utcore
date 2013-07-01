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
 * Class to bind a function to a parameter.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_DETAIL_BINDER_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_DETAIL_BINDER_H_INCLUDED__
 
#include "ResultVector.h" 
#include "ResultMatrix.h" 
 
namespace Ubitrack { namespace Math { namespace Function { namespace Detail {

/**
 * \internal
 * Binds a function to a parameter.
 * Multiple parameters are chained in this fashion:
 * <pre>Binder< Binder< Func, Param1 >, Param2 ></pre>
 */
template< class CFunc, class CParam >
class Binder
{
public:
	Binder( const CFunc& f, const CParam& p )
		: m_func( f )
		, m_param( p )
		, m_result( f.size() )
	{}

	std::size_t size() const
	{ return m_func.size(); }

private:
	template< class, class > friend class Binder;
	static const std::size_t staticSize = CFunc::staticSize;
	static const bool wantsJacobian = CFunc::wantsJacobian || CParam::wantsJacobian;
	
	template< class ParameterVector >
	const ResultVector< CFunc::staticSize >& value( const ParameterVector& ) const
	{ return m_result; }


	// internal evaluation functions
	
	// one-parameter function
	template< class ParameterVector, class DestinationVector >
	void i_evaluate( const ParameterVector& p, DestinationVector& d ) const
	{
		m_param.i_evaluateInternal( p );
		m_func.i_evaluate( p, d, m_param.value( p ) );
	}
	
	// two-parameter function
	template< class ParameterVector, class DestinationVector, class Param1 >
	void i_evaluate( const ParameterVector& p, DestinationVector& d, const Param1& p1 ) const
	{
		m_param.i_evaluateInternal( p );
		m_func.i_evaluate( p, d, m_param.value( p ), p1 );
	}
	
	// three-parameter function
	template< class ParameterVector, class DestinationVector, class Param1, class Param2 >
	void i_evaluate( const ParameterVector& p, DestinationVector& d, const Param1& p1, const Param2& p2 ) const
	{
		m_param.i_evaluateInternal( p );
		m_func.i_evaluate( p, d, m_param.value( p ), p1, p2 );
	}

	// only evaluate parameters
	template< class ParameterVector >
	void i_evaluateParameters( const ParameterVector& p ) const
	{
		m_param.i_evaluateInternal( p );
		m_func.i_evaluateParameters( p );
	}

	
	// internal evaluation with internal storage
	
	// evaluate and store in internal result vector
	template< class ParameterVector >
	void i_evaluateInternal( const ParameterVector& p ) const
	{
		i_evaluate( p, m_result );
	}

	
	// internal jacobian function
	
	// multiply jacobian of current parameter (no right-hand-side parameters)
	template< std::size_t LHSize, class ParameterVector, class LeftHand, class DestinationMatrix >
	void i_multiplyJacobian( const ParameterVector& p, const LeftHand& l, DestinationMatrix& j ) const
	{
		// assume evaluate() has already been called
		// first, pass on to jacobian of parameter to left
#ifdef __clang__
		m_func.template i_multiplyJacobian< LHSize >( p, l, j, m_param.value( p ) );
#else
		m_func.i_multiplyJacobian< LHSize >( p, l, j, m_param.value( p ) );
#endif
		if ( CParam::wantsJacobian )
		{
			Detail::ResultMatrix< LHSize * CParam::staticSize, LHSize, CParam::staticSize > jResult( j.size1(), m_param.size() );
			m_func.i_multiplyJacobian1( p, l, jResult, m_param.value( p ) );
#ifdef __clang__
			m_param.template i_multiplyJacobian< LHSize >( p, jResult, j );
#else
			m_param.i_multiplyJacobian< LHSize >( p, jResult, j );
#endif
		}
	}
	
	// multiply jacobian of current parameter (one right-hand-side parameter)
	template< std::size_t LHSize, class ParameterVector, class LeftHand, class DestinationMatrix, class Param1 >
	void i_multiplyJacobian( const ParameterVector& p, const LeftHand& l, DestinationMatrix& j, const Param1& p1 ) const
	{
		// assume evaluate() has already been called
		// first, pass on to jacobian of parameter to left
#ifdef __clang__
		m_func.template i_multiplyJacobian< LHSize >( p, l, j, m_param.value( p ), p1 );
#else
		m_func.i_multiplyJacobian< LHSize >( p, l, j, m_param.value( p ), p1 );
#endif 
		
		if ( CParam::wantsJacobian )
		{
			Detail::ResultMatrix< LHSize * CParam::staticSize, LHSize, CParam::staticSize > jResult( j.size1(), m_param.size() );
			m_func.i_multiplyJacobian1( p, l, jResult, m_param.value( p ), p1 );
#ifdef __clang__
			m_param.template i_multiplyJacobian< LHSize >( p, jResult, j );
#else
			m_param.i_multiplyJacobian< LHSize >( p, jResult, j );
#endif
		}
	}
	
	// multiply jacobian of current parameter (two right-hand-side parameters)
	template< std::size_t LHSize, class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void i_multiplyJacobian( const ParameterVector& p, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2 ) const
	{
		// assume evaluate() has already been called
		// first, pass on to jacobian of parameter to left
		// m_func.i_multiplyJacobian< LHSize >( p, l, j, m_param.value( p ), p1, p2 ); // only 3 parameters supported
		
		if ( CParam::wantsJacobian )
		{
			Detail::ResultMatrix< LHSize * CParam::staticSize, LHSize, CParam::staticSize > jResult( j.size1(), m_param.size() );
			m_func.i_multiplyJacobian1( p, l, jResult, m_param.value( p ), p1, p2 );
			m_param.i_multiplyJacobian< LHSize >( p, jResult, j );
		}
	}

	// pass-through jacobian multiplication (one right-hand-side parameter)
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1 >
	void i_multiplyJacobian1( const ParameterVector& p, const LeftHand& l, DestinationMatrix& j, const Param1& p1 ) const
	{
		m_func.i_multiplyJacobian2( p, l, j, m_param.value( p ), p1 );
	}

	// pass-through jacobian multiplication (1st of two right-hand-side parameter)
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void i_multiplyJacobian1( const ParameterVector& p, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2 ) const
	{
		m_func.i_multiplyJacobian2( p, l, j, m_param.value( p ), p1, p2 );
	}

	// pass-through jacobian multiplication (2nd of two right-hand-side parameter)
	template< class ParameterVector, class LeftHand, class DestinationMatrix, class Param1, class Param2 >
	void i_multiplyJacobian2( const ParameterVector& p, const LeftHand& l, DestinationMatrix& j, const Param1& p1, const Param2& p2 ) const
	{
		m_func.i_multiplyJacobian3( p, l, j, m_param.value( p ), p1, p2 );
	}


public:	
	// outermost evaluation functions
	
	// function evaluation
	template< class ParameterVector, class DestinationVector >
	void evaluate( const ParameterVector& p, DestinationVector& d ) const
	{
		i_evaluate( p, d );
	}

	// top-level jacobian call
	template< class ParameterVector, class DestinationMatrix >
	void jacobian( const ParameterVector& p, DestinationMatrix& j ) const
	{
		assert( staticSize == 0 || j.size1() == staticSize );
		i_evaluateParameters( p ); 
		i_multiplyJacobian< staticSize >( p, ublas::identity_matrix< double >( size() ), j );
	}
	
	// top-level jacobian call
	template< class ParameterVector, class ResultVector, class DestinationMatrix >
	void evaluateWithJacobian( const ParameterVector& p, ResultVector& r, DestinationMatrix& j ) const
	{
		assert( staticSize == 0 || j.size1() == staticSize );
		i_evaluate( p, r ); 
		i_multiplyJacobian< staticSize >( p, ublas::identity_matrix< double >( size() ), j );
	}
	
	
private:
	CFunc m_func;
	CParam m_param;
	mutable ResultVector< CFunc::staticSize > m_result;
};


} } } } // namespace Ubitrack::Math::Function::Detail

#endif

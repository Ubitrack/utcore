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
 * Main include file of the new multivariate function framework.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_FUNCTION_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_FUNCTION_H_INCLUDED__

#include "Detail/Parameter.h"
#include "Detail/FixedParameterRef.h"
#include "Detail/FixedParameterCopy.h"
#include "Detail/Binder.h"
 
namespace Ubitrack { namespace Math { namespace Optimization { namespace Function {

namespace Detail
{
	/** 
	 * \internal
	 * used internally to distinguish parameter classes from other classes
	 */
	template< class P >
	class ParameterWrapper
		: public P
	{
	public:
		ParameterWrapper( const P& _p )
			: P( _p )
		{}
	};

} // namespace Detail


/** creates a parameter object to refer to a parameter that should be optimized */
template< std::size_t Size >
Detail::ParameterWrapper< Detail::Parameter< Size > > parameter( std::size_t nStart )
{ return Detail::ParameterWrapper< Detail::Parameter< Size > >( Detail::Parameter< Size >( nStart ) ); }

/** 
 * creates a parameter object to refer to a constant (non-optimized) parameter. 
 * The parameter value reference must stay constant during the lifetime of the object.
 */
template< std::size_t Size, class CVector >
Detail::ParameterWrapper< Detail::FixedParameterRef< Size, CVector > > fixedParameterRef( const CVector& v )
{ return Detail::ParameterWrapper< Detail::FixedParameterRef< Size, CVector > >( Detail::FixedParameterRef< Size, CVector >( v ) ); }

/** 
 * creates a parameter object to refer to a constant (non-optimized) parameter. 
 * The parameter value is copied to an internal vector.
 */
template< std::size_t Size, class CVector >
Detail::ParameterWrapper< Detail::FixedParameterCopy< Size > > fixedParameterCopy( const CVector& v )
{ return Detail::ParameterWrapper< Detail::FixedParameterCopy< Size > >( Detail::FixedParameterCopy< Size >( v ) ); }


/**
 * Bind a function object to a (final) parameter
 */
template< class CFunc, class CParam >
Detail::Binder< CFunc, CParam > operator<<( const CFunc& func, const Detail::ParameterWrapper< CParam >& param )
{ return Detail::Binder< CFunc, CParam >( func, param ); }

/**
 * Bind a function object to a parameter computed by another function
 */
template< class CFunc, class CBFunc, class CBParam >
Detail::Binder< CFunc, Detail::Binder< CBFunc, CBParam > > operator<<( const CFunc& func, const Detail::Binder< CBFunc, CBParam >& binder )
{ return Detail::Binder< CFunc, Detail::Binder< CBFunc, CBParam > >( func, binder ); }


}}}} // namespace Ubitrack::Math::Optimization::Function

#endif

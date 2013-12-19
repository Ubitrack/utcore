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
 * class for parameters which are optimized (i.e. part of the parameter vector)
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_MATH_FUNCTION_DETAIL_PARAMETER_H_INCLUDED__
#define __UBITRACK_MATH_FUNCTION_DETAIL_PARAMETER_H_INCLUDED__


#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace Ubitrack { namespace Math { namespace Optimization { namespace Function { namespace Detail {

/**
 * class for parameters which are optimized (i.e. part of the parameter vector)
 */
template< std::size_t Size = 0 >
class Parameter
{
public:
	Parameter( std::size_t nStart, std::size_t nSize = Size )
		: m_range( nStart, nStart + nSize )
	{
		assert( !Size || nSize == Size );
	}

	std::size_t size() const
	{ return m_range.size(); }

private:	
	template< class, class > friend class Binder;
	
	static const std::size_t staticSize = Size;
	static const bool wantsJacobian = true;

	template< class ParameterVector >
	const boost::numeric::ublas::vector_range< const ParameterVector > value( const ParameterVector& p ) const
	{ return boost::numeric::ublas::vector_range< const ParameterVector >( p, m_range ); }

	/** store jacobian in destination matrix */
	template< std::size_t LHSize, class ParameterVector, class LeftHand, class DestinationMatrix >
	void i_multiplyJacobian( const ParameterVector&, const LeftHand& l, DestinationMatrix& j ) const
	{
		boost::numeric::ublas::matrix_range< DestinationMatrix > dest( j, boost::numeric::ublas::range( 0, j.size1() ), m_range );
		dest = l;
	}

	template< class ParameterVector >
	void i_evaluateInternal( const ParameterVector& ) const
	{}
	
	boost::numeric::ublas::range m_range;
};

}}}}} // namespace Ubitrack::Math::Optimization::Function::Detail

#endif

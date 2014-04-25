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
 * @ingroup tracking_algorithms
 * @file
 * Functions for non-linearly optimized tooltip/hotspot calibration.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

#ifndef __UBITRACK_ALGROITHM_TOOLTIP_OPTIMIZATION_H_INCLUDED__
#define __UBITRACK_ALGROITHM_TOOLTIP_OPTIMIZATION_H_INCLUDED__

#include "TipCalibration.h" // includes std::vector/Pose/Vector
#include <utMath/Blas1.h>
#include <utMath/Optimization/Optimization.h>
#include <utMath/Optimization/LevenbergMarquardt.h>


// Boost
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

namespace Ubitrack { namespace Algorithm { namespace ToolTip {

/**
 * @internal 
 * A minimization function for the non-linear optimization
 * can be applied with Levenberg-Marquardt
 */
template< typename InputIterator >
class MultiplePoseSinglePointTransformation
{

	typedef typename std::iterator_traits< InputIterator >::value_type::value_type VType;
	
protected:
	const InputIterator m_iBegin;
	const InputIterator m_iEnd;
	
public:
	/** 
	 * constructor.
	 * @param iBegin iterator to the beginning of a container with projections(must stay constant during lifetime of the object)
	 * @param iEnd iterator to the end of a container with projections(must stay constant during lifetime of the object)
	 */
	MultiplePoseSinglePointTransformation( const InputIterator iBegin, const InputIterator iEnd )
		: m_iBegin( iBegin )
		, m_iEnd( iEnd )
	{}

	/**
	 * return the size of the result vector
	 */
	std::size_t size() const
	{ return ( std::distance( m_iBegin, m_iEnd ) ); }

	/**
	 * @param result 6-vector to store the result in
	 * @param input containing the tooltip parameters (tooltip in worldand tooltip offset)
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		std::size_t i ( 0 );
		for ( InputIterator it ( m_iBegin ); it != m_iEnd; ++i, ++it )
		{
			const Ubitrack::Math::Vector< VType, 3 > vec1( input[ 0 ], input[ 1 ], input[ 2 ] );
			const Ubitrack::Math::Vector< VType, 3 > vec2( input[ 3 ], input[ 4 ], input[ 5 ] );
			const Math::Vector< VType, 3 > vec = vec1 - ( (*it) * vec2 );
			result[ i ] = Ubitrack::Math::norm_2( vec );
		}
		// std::cout << "Result : \n" << result << std::endl;
	}
	
	/**
	 * @param result vector to store the result in
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the Jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		// TODO: implement as one function (more efficient)
		evaluate( result, input );
		jacobian( input, J );
	}

	/**
	 * @param input containing the parameters (to be optimized)
	 * @param J matrix to store the Jacobian (evaluated for input) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		
		const VType x1 = input[ 0 ];
		const VType y1 = input[ 1 ];
		const VType z1 = input[ 2 ];
		const VType x2 = input[ 3 ];
		const VType y2 = input[ 4 ];
		const VType z2 = input[ 5 ];
		std::size_t i( 0 );
		for ( InputIterator it ( m_iBegin ); it != m_iEnd; ++i, ++it )
		{
			const VType qx = it->rotation().x();
			const VType qy = it->rotation().y();
			const VType qz = it->rotation().z();
			const VType qw = it->rotation().w();
			const VType tx = it->translation()[ 0 ];
			const VType ty = it->translation()[ 1 ];
			const VType tz = it->translation()[ 2 ];
			
			const VType t2 = qy*qy;
			const VType t3 = t2*2.0;
			const VType t4 = qz*qz;
			const VType t5 = t4*2.0;
			const VType t6 = t3+t5-1.0;
			const VType t7 = t6*x2;
			const VType t8 = qw*qz*2.0;
			const VType t15 = qx*qy*2.0;
			const VType t9 = t8-t15;
			const VType t10 = t9*y2;
			const VType t11 = qw*qy*2.0;
			const VType t12 = qx*qz*2.0;
			const VType t13 = t11+t12;
			const VType t16 = t13*z2;
			const VType t14 = t7+t10-t16-tx+x1;
			const VType t17 = fabs(t14);
			const VType t19 = qx*qx;
			const VType t20 = t19*2.0;
			const VType t21 = qw*qx*2.0;
			const VType t22 = qy*qz*2.0;
			const VType t24 = t5+t20-1.0;
			const VType t25 = t24*y2;
			const VType t26 = t8+t15;
			const VType t27 = t26*x2;
			const VType t28 = t21-t22;
			const VType t29 = t28*z2;
			const VType t30 = t25-t27+t29-ty+y1;
			const VType t18 = fabs(t30);
			const VType t33 = t3+t20-1.0;
			const VType t34 = t33*z2;
			const VType t35 = t11-t12;
			const VType t36 = t35*x2;
			const VType t37 = t21+t22;
			const VType t38 = t37*y2;
			const VType t39 = t34+t36-t38-tz+z1;
			const VType t23 = fabs(t39);
			const VType t31 = t17*t17;
			const VType t32 = t18*t18;
			const VType t40 = t23*t23;
			const VType t41 = t31+t32+t40;
			const VType t42 = 1.0/sqrt(t41);
			const VType t43 = (t14/fabs(t14));
			const VType t44 = (t30/fabs(t30));
			const VType t45 = (t39/fabs(t39));
			J( i, 0 ) = t17*t42*t43;
			J( i, 1 ) = t18*t42*t44;
			J( i, 2 ) = t23*t42*t45;
			J( i, 3 ) = t42*(t6*t17*t43*2.0-t18*t26*t44*2.0+t23*t35*t45*2.0)*(1.0/2.0);
			J( i, 4 ) = t42*(t9*t17*t43*2.0+t18*t24*t44*2.0-t23*t37*t45*2.0)*(1.0/2.0);
			J( i, 5 ) = t42*(t13*t17*t43*-2.0+t18*t28*t44*2.0+t23*t33*t45*2.0)*(1.0/2.0);
		}
		// std::cout << "input  :\n" << input << "\n";
		// std::cout << "J      :\n" << J.size1() << " x " << J.size2() << std::endl;
	}
};



bool estimatePosition3D_6D( Math::Vector< double, 3 >& pw
	, const std::vector< Math::Pose >& poses 
	, Math::Vector< double, 3 >& pm
	, const Math::Optimization::OptTerminate& criteria )
	{
	
		const std::size_t n ( poses.size() );
		
		// 1) estimate a first initial guess
		if( !estimatePosition3D_6D( pw, poses, pm ) )
			return false;
		
		// 2) set some initial values to be optimized non-linearly
		Math::Vector< double, 6 > paramVector;
		paramVector( 0 ) = pw[ 0 ];
		paramVector( 1 ) = pw[ 1 ];
		paramVector( 2 ) = pw[ 2 ];
		paramVector( 3 ) = pm[ 0 ];
		paramVector( 4 ) = pm[ 1 ];
		paramVector( 5 ) = pm[ 2 ];

		
		// 3) prepare the expectation values of the minimization function
		Math::Vector< double > measurement = Math::Vector< double >::zeros( n );
		
		// 4) set the evaluation function
		ToolTip::MultiplePoseSinglePointTransformation< std::vector< Math::Pose >::const_iterator > func( poses.begin(), poses.end() );
		
		// 5) perform optimization
		double residual = Ubitrack::Math::Optimization::levenbergMarquardt( func, paramVector, measurement, criteria, Math::Optimization::OptNoNormalize() );	
		// std::cout << "Residual Error " << residual << std::endl;
		pw[ 0 ] = paramVector[ 0 ];
		pw[ 1 ] = paramVector[ 1 ] ;
		pw[ 2 ] = paramVector[ 2 ] ;
		pm[ 0 ] = paramVector[ 3 ] ;
		pm[ 1 ] = paramVector[ 4 ] ;
		pm[ 2 ] = paramVector[ 5 ] ;
		
		return true;
	}


}}} // namespace Ubitrack::Algorithm::ToolTip

#endif // __UBITRACK_ALGROITHM_TOOLTIP_OPTIMIZATION_H_INCLUDED__
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
 *
 * Math::Vector< 10 > wrapper to explicitly represent a quadric
 * (or quadratic surface) of the following form:
 * @f$ax^2+by^2+cz^2+2fyz+2gzx+2hxy+2px+2qy+2rz+d=0@f$.
 *
 * The header includes functors for common operations on quadrics.
 * The functors can easily be applied to STL-containers like 
 * std::vectors, std::list, etc, that contain 10-vectors 
 * representing quadratic equations or explicit quadrics.
 *
 * For some general information, please have a look at:
 * http://mathworld.wolfram.com/QuadraticSurface.html
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __H__QUADRIC_FUNCTORS__
#define __H__QUADRIC_FUNCTORS__

// Ubitrack
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utUtil/Exception.h>
#include <utMath/Functors/MatrixFunctors.h>

// std
#include <numeric>
#include <algorithm>
#include <functional>

namespace Ubitrack { namespace Math { namespace Geometry {

/**
 * @ingroup math
 * Wraps a Math::Vector< 10, T > to explicitly represent a quadric.
 *
 * @tparam T type of conic ( e.g. \c double or \c float )
 */
template< typename T >
struct Quadric
	: public Math::Vector< 10, T >
{
public:
	/** Default Constructor */
	Quadric( ){}

	/** Copy Constructor for implicit conversion */
	Quadric( const Math::Vector< 10, T >& quadric )
		: Math::Vector< 10, T >()
	{
		(*this)[ 0 ] = quadric[ 0 ];// a (1st semi-axis for ellipsoids)
		(*this)[ 1 ] = quadric[ 1 ];// b (2nd semi-axis for ellipsoids)
		(*this)[ 2 ] = quadric[ 2 ];// c (3rd semi-axis for ellipsoids)
		(*this)[ 3 ] = quadric[ 3 ];// f (= 0 for ellipsoids)
		(*this)[ 4 ] = quadric[ 4 ];// g (= 0 for ellipsoids)
		(*this)[ 5 ] = quadric[ 5 ];// h (= 0 for ellipsoids)
		(*this)[ 6 ] = quadric[ 6 ];// p (= x for ellipsoids)
		(*this)[ 7 ] = quadric[ 7 ];// q (= y for ellipsoids)
		(*this)[ 8 ] = quadric[ 8 ];// r (= z for ellipsoids)
		(*this)[ 9 ] = quadric[ 9 ];// d (-1/0/1, -1 for ellipsoids )
	}
};


/**
 * @ingroup math
 * Projects a quadric onto the image plane by a 3x4 projection matrix.
 *
 * This functor class can be applied to STL-containers
 * of quadrics via a STL-algorithms.
 */
template< typename T >
struct ProjectQuadric
	: public std::binary_function< Math::Matrix< 3, 4, T >, Math::Vector< 10, T >, Math::Vector< 6, T > >
{
public:
	/**
	 * @ingroup math
	 * Projects a quadric resulting in a conic.
	 *
	 * Usually one can image a linear system of \f$ (P * Q * P^T)^{-1} \f$ with P as
	 * a 3x4 projection matrix and Q the symmetric 4x4 matrix constructed from the
	 * 10-vector representation.
	 *
	 * @tparam T type of quadric ( e.g. \c double or \c float )
	 * @param projection the 3x4 projection matrix 
	 * @param quadric a 10-vector including the quadric's parameters
	 * @return the resulting conic
	 */
	Math::Vector< 6, T > operator() ( const Math::Matrix< 3, 4, T > &projection, const Math::Vector< 10, T > &quadric ) const
	{
		const T a ( quadric( 0 ) );
		const T b ( quadric( 1 ) );
		const T c ( quadric( 2 ) );
		const T d ( quadric( 9 ) );
	
		const T f ( quadric( 3 ) );
		const T g ( quadric( 4 ) );
		const T h ( quadric( 5 ) );
	
		const T p ( quadric( 6 ) ); // or x
		const T q ( quadric( 7 ) ); // or y
		const T r ( quadric( 8 ) ); // or z
	
		const T p11 ( projection( 0, 0 ) );
		const T p12 ( projection( 0, 1 ) );
		const T p13 ( projection( 0, 2 ) );
		const T p14 ( projection( 0, 3 ) );
		const T p21 ( projection( 1, 0 ) );
		const T p22 ( projection( 1, 1 ) );
		const T p23 ( projection( 1, 2 ) );
		const T p24 ( projection( 1, 3 ) );
		const T p31 ( projection( 2, 0 ) );
		const T p32 ( projection( 2, 1 ) );
		const T p33 ( projection( 2, 2 ) );
		const T p34 ( projection( 2, 3 ) );

		const T conic_a = p11*(p11*a+f*p12+g*p13+p*p14)+p12*(f*p11+p12*b+h*p13+q*p14)+p13*(g*p11+h*p12+p13*c+r*p14)+p14*(p11*p+p12*q+p13*r+d*p14);
		const T conic_b = 2*(p11*(p21*a+f*p22+g*p23+p*p24)+p12*(f*p21+p22*b+h*p23+q*p24)+p13*(g*p21+h*p22+p23*c+r*p24)+p14*(p21*p+p22*q+p23*r+d*p24));
		const T conic_c = p21*(p21*a+f*p22+g*p23+p*p24)+p22*(f*p21+p22*b+h*p23+q*p24)+p23*(g*p21+h*p22+p23*c+r*p24)+p24*(p21*p+p22*q+p23*r+d*p24);
		const T conic_d = 2*(p11*(p31*a+f*p32+g*p33+p*p34)+p12*(f*p31+p32*b+h*p33+q*p34)+p13*(g*p31+h*p32+p33*c+r*p34)+p14*(p31*p+p32*q+p33*r+d*p34));
		const T conic_e = 2*(p21*(p31*a+f*p32+g*p33+p*p34)+p22*(f*p31+p32*b+h*p33+q*p34)+p23*(g*p31+h*p32+p33*c+r*p34)+p24*(p31*p+p32*q+p33*r+d*p34));
		const T conic_f = p31*(p31*a+f*p32+g*p33+p*p34)+p32*(f*p31+p32*b+h*p33+q*p34)+p33*(g*p31+h*p32+p33*c+r*p34)+p34*(p31*p+p32*q+p33*r+d*p34);

		// make it a point-conic:
		const T divisor =  1/( (conic_a*(conic_e*conic_e)+conic_c*(conic_d*conic_d)+(conic_b*conic_b)*conic_f-conic_a*conic_c*conic_f*static_cast< T >( 4 )-conic_b*conic_d*conic_e) );
		
		Math::Vector< 6, T > iConic;
		iConic( 0 ) = -(conic_c*conic_f*static_cast< T >( 4 ) -conic_e*conic_e) * divisor;
		iConic( 1 ) =  2*(conic_b*conic_f*static_cast< T >( 2 ) -conic_d*conic_e) * divisor;
		iConic( 2 )  = -(conic_a*conic_f*static_cast< T >( 4 ) -conic_d*conic_d) * divisor;
		iConic( 3 )  = 2*(-(conic_b*conic_e-conic_c*conic_d*static_cast< T >( 2 ) ) * divisor );
		iConic( 4 )  = 2*(conic_a*conic_e*static_cast< T >( 2 ) -conic_b*conic_d) * divisor;
		iConic( 5 )  = -(conic_a*conic_c*static_cast< T >( 4 ) -conic_b*conic_b) * divisor;	
		// or just make it a line-conic:
		/*
		Math::Vector< 6, T > iConic;
		iConic( 0 ) = conic_a;
		iConic( 1 ) = conic_b;
		iConic( 2 ) = conic_c;
		iConic( 3 ) = conic_d;
		iConic( 4 ) = conic_e;
		iConic( 5 ) = conic_f;
		*/
		return Math::Vector< 6, T >( iConic );
	}
};

/**
 * @ingroup math
 * Generates a quadric in general description(10-vector) from a 6-vector.
 *
 * This functor class can be applied to STL-containers
 * of ellipsoids via a STL-algorithms.
 */
template< typename T >
struct Ellipsoid2Quadric
	: public std::unary_function< Math::Vector< 6, T >, Math::Vector< 10, T > >
{
public:
	/**
	 * @ingroup math
	 * Projects a quadric resulting in a conic.
	 *
	 * @tparam T type of quadric ( e.g. \c double or \c float )
	 * @param projection the 3x4 projection matrix 
	 * @param quadric a 10-vector including the parameters of the quadric
	 * @return the resulting conic
	 */
	Math::Vector< 10, T > operator() ( const Math::Vector< 6, T > &ellipsoid ) const
	{
		const T a ( ellipsoid( 0 ) ); //1st semi-axis
		const T b ( ellipsoid( 1 ) ); //2nd semi-axis
		const T c ( ellipsoid( 2 ) ); //3rd semi-axis
		const T d ( 1 ); //just fo nicer code writing ;)
	
		const T p ( ellipsoid( 3 ) ); // or x-position
		const T q ( ellipsoid( 4 ) ); // or y-position
		const T r ( ellipsoid( 5 ) ); // or z-position
		// we assume H * D * H'
		// D := diagonal matrix with d1, d2, d3, d4
		// for ellipsoids d1 = d2 = d3 = 1 and d4 = -1
		// H := quadric matrix with rows: [a f g p], [0 b h q], [0 0 c r], [0 0 0 1]
		// for ellipsoids f = g = h = 0;
		// results in symmetric matrix:
		//[ d1*a^2 + d3*g^2 + d2*h^2 + d4*p^2, b*d2*h + d3*f*g + d4*p*q, c*d3*g + d4*p*r, d*d4*p]
		//[          b*d2*h + d3*f*g + d4*p*q, d2*b^2 + d3*f^2 + d4*q^2, c*d3*f + d4*q*r, d*d4*q]
		//[                   c*d3*g + d4*p*r,          c*d3*f + d4*q*r, d3*c^2 + d4*r^2, d*d4*r]
		//[                            d*d4*p,                   d*d4*q,          d*d4*r, d^2*d4]

		Math::Vector< 10, T > quadric;
		// diagonal entries first
		quadric( 0 ) = (a*a)-(p*p);
		quadric( 1 ) = (b*b)-(q*q);
		quadric( 2 ) = (c*c)-(r*r);
		quadric( 9 ) = -(d*d);
		// symmetric part
		quadric( 3 ) = -(p*q);
		quadric( 4 ) = -(p*r);
		quadric( 5 ) = -(q*r);
		// last column/lowest row
		quadric( 6 ) = -(d*p);
		quadric( 7 ) = -(d*q);
		quadric( 8 ) = -(d*r);

		return Math::Vector< 10, T >( quadric );
	}
};


/**
 * @ingroup math
 * Projects an ellipsoid onto the image plane by a 3x4 projection matrix.
 *
 * This functor class can be applied to STL-containers
 * of ellipsoids via a STL-algorithms.
 */
template< typename T >
struct ProjectEllipsoid
	: public std::binary_function< Math::Matrix< 3, 4, T >, Math::Vector< 6, T >, Math::Vector< 6, T > >
{
public:
	/**
	 * @ingroup math
	 * Projects an ellipsoid onto the image plane by a 3x4 projection matrix.
	 *
	 * @tparam T type of ellipsoid ( e.g. \c double or \c float )
	 * @param projection the 3x4 projection matrix 
	 * @param ellipsoid a 6-vector including the parameters of the ellipsoid
	 * @return the resulting conic
	 */
	Math::Vector< 6, T > operator() ( const Math::Matrix< 3, 4, T > &projection, const Math::Vector< 6, T > &ellipsoid ) const
	{
		const T a ( ellipsoid( 0 ) ); //1st semi-axis
		const T b ( ellipsoid( 1 ) ); //2nd semi-axis
		const T c ( ellipsoid( 2 ) ); //3rd semi-axis
		//const T d ( 1 ); //not necessary
	
		const T p ( ellipsoid( 3 ) ); // or x
		const T q ( ellipsoid( 4 ) ); // or y
		const T r ( ellipsoid( 5 ) ); // or z
	
		const T p11 ( projection( 0, 0 ) );
		const T p12 ( projection( 0, 1 ) );
		const T p13 ( projection( 0, 2 ) );
		const T p14 ( projection( 0, 3 ) );
		const T p21 ( projection( 1, 0 ) );
		const T p22 ( projection( 1, 1 ) );
		const T p23 ( projection( 1, 2 ) );
		const T p24 ( projection( 1, 3 ) );
		const T p31 ( projection( 2, 0 ) );
		const T p32 ( projection( 2, 1 ) );
		const T p33 ( projection( 2, 2 ) );
		const T p34 ( projection( 2, 3 ) );
	
	
		const T conic_a = p11*p11*a*a+p12*p12*b*b+p13*p13*c*c+(-p11*p-p12*q-p13*r-p14)*(p11*p+p12*q+p13*r+p14);
		const T conic_b = 2*(p11*a*a*p21+p12*b*b*p22+p13*c*c*p23+(-p11*p-p12*q-p13*r-p14)*(p21*p+p22*q+p23*r+p24));
		const T conic_c = p21*p21*a*a+p22*p22*b*b+p23*p23*c*c+(-p21*p-p22*q-p23*r-p24)*(p21*p+p22*q+p23*r+p24);
		const T conic_d = 2*(p11*a*a*p31+p12*b*b*p32+p13*c*c*p33+(-p11*p-p12*q-p13*r-p14)*(p31*p+p32*q+p33*r+p34));
		const T conic_e = 2*(p21*a*a*p31+p22*b*b*p32+p23*c*c*p33+(-p21*p-p22*q-p23*r-p24)*(p31*p+p32*q+p33*r+p34));
		const T conic_f = p31*p31*a*a+p32*p32*b*b+p33*p33*c*c+(-p31*p-p32*q-p33*r-p34)*(p31*p+p32*q+p33*r+p34);

		
		// make it a point-conic:
		const T divisor =  1/( (conic_a*(conic_e*conic_e)+conic_c*(conic_d*conic_d)+(conic_b*conic_b)*conic_f-conic_a*conic_c*conic_f*static_cast< T >( 4 )-conic_b*conic_d*conic_e) );

		Math::Vector< 6, T > iConic;		
		iConic( 0 ) = -(conic_c*conic_f*static_cast< T >( 4 ) -conic_e*conic_e) * divisor;
		iConic( 1 ) =  2*(conic_b*conic_f*static_cast< T >( 2 ) -conic_d*conic_e) * divisor;
		iConic( 2 )  = -(conic_a*conic_f*static_cast< T >( 4 ) -conic_d*conic_d) * divisor;
		iConic( 3 )  = 2*(-(conic_b*conic_e-conic_c*conic_d*static_cast< T >( 2 ) ) * divisor );
		iConic( 4 )  = 2*(conic_a*conic_e*static_cast< T >( 2 ) -conic_b*conic_d) * divisor;
		iConic( 5 )  = -(conic_a*conic_c*static_cast< T >( 4 ) -conic_b*conic_b) * divisor;	
		
		// or just make it a line-conic:
		/*
		Math::Vector< 6, T > iConic;
		iConic[ 0 ] = conic_a;
		iConic[ 1 ] = conic_b;
		iConic[ 2 ] = conic_c;
		iConic[ 3 ] = conic_d;
		iConic[ 4 ] = conic_e;
		iConic[ 5 ] = conic_f;
		*/
		return Math::Vector< 6, T >( iConic );
	}
};

/**
 * @ingroup math
 * Projects a spheroid onto the image plane by a 3x4 projection matrix.
 *
 * This functor class can be applied to STL-containers
 * of spheroids via a STL-algorithms.
 */
template< typename T >
struct ProjectSpheroid
	: public std::binary_function< Math::Matrix< 3, 4, T >, Math::Vector< 4, T >, Math::Vector< 6, T > >
{
protected:
	const ProjectEllipsoid< T > m_Projector;

public:

	/** Standard constructor */
	ProjectSpheroid()
		: std::binary_function< Math::Matrix< 3, 4, T >, Math::Vector< 4, T >, Math::Vector< 6, T > >()
		, m_Projector()
		{}

	/**
	 * @ingroup math
	 * Projects an spheroid onto the image plane by a 3x4 projection matrix.
	 *
	 * @tparam T type of spheroid ( e.g. \c double or \c float )
	 * @param projection the 3x4 projection matrix 
	 * @param ellipsoid a 6-vector including the parameters of the ellipsoid
	 * @return the resulting conic
	 */
	Math::Vector< 6, T > operator() ( const Math::Matrix< 3, 4, T > &projection, const Math::Vector< 4, T > &spheroid ) const
	{
		Math::Vector< 6, T > ellipsoid;
		ellipsoid( 0 ) = spheroid( 0 ); // 1st semi axis
		ellipsoid( 1 ) = spheroid( 0 ); // 2nd semi axis
		ellipsoid( 2 ) = spheroid( 1 ); // 3rd semi axis
		ellipsoid( 3 ) = spheroid( 2 ); // x-position
		ellipsoid( 4 ) = spheroid( 3 ); // y-position
		ellipsoid( 5 ) = spheroid( 1 ); // z-position (similar to 3rd semi axis)
		return m_Projector( projection, ellipsoid );
	}
};


}}} // namespace Ubitrack::Math::Geometry


#endif // __H__QUADRIC__

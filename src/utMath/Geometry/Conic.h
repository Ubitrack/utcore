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
 * Math::Vector< 6 > wrapper to explicitly represent conics.
 *
 * The header includes functors for common operations on conics.
 * The functors can easily be applied to STL-containers like 
 * std::vectors, std::list, etc, that contain 6-vectors 
 * representing conics or explicit conics.
 *
 * For some general information, please have a look at:
 * http://mathworld.wolfram.com/QuadraticCurve.html
 * http://mathworld.wolfram.com/ConicSection.html
 * http://mathworld.wolfram.com/Ellipse.html
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 


#ifndef __H__CONIC_FUNCTORS__
#define __H__CONIC_FUNCTORS__

// Ubitrack
#include <utMath/Vector.h>
#include <utMath/Matrix.h>
#include <utUtil/Exception.h>
#include <utUtil/StaticAssert.h>
#include <utMath/Functors/MatrixFunctors.h>

// std
#include <numeric>
#include <algorithm>
#include <functional>

#ifndef M_PI
#define _USE_MATH_DEFINES //for having PI
#include <math.h>
#endif

namespace Ubitrack { namespace Math { namespace Geometry {


/**
 * @ingroup math
 * Wraps a Math::Vector< T, 6 > to explicitly represent a conic.
 *
 * @tparam T type of conic ( e.g. \c double or \c float )
 */
template< typename T >
struct Conic
	: public Math::Vector< T, 6 >
{
public:
	/** Default Constructor */
	Conic( ){}

	/** Copy Constructor for implicit conversion */
	Conic( const Math::Vector< T, 6 >& conic )
		: Math::Vector< T, 6 >()
	{
		(*this)[ 0 ] = conic[ 0 ];
		(*this)[ 1 ] = conic[ 1 ];
		(*this)[ 2 ] = conic[ 2 ];
		(*this)[ 3 ] = conic[ 3 ];
		(*this)[ 4 ] = conic[ 4 ];
		(*this)[ 5 ] = conic[ 5 ];
	}
};


/**
 * @ingroup math
 * Changes the representation of a conic from vectorial to matrix.
 * 
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct MatrixFromConic
	: public std::unary_function< Math::Vector< T, 6 >, Math::Matrix< T, 3, 3 > >
{
public:
	/**
	 * @ingroup math
	 * Changes the representation of a conic from vectorial to a matrix.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the conic expressed as a vector
	 * @return the 3x3-matrix representation of the conic
	 */
	Math::Matrix< T, 3, 3 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major > matrix( 3, 3 );
		matrix( 0, 0 ) = conic[ 0 ];
		matrix( 1, 0 ) = matrix( 0, 1 ) = conic[ 1 ] * static_cast< T > ( 0.5 );
		matrix( 1, 1 ) = conic[ 2 ];
		matrix( 2, 0 ) = matrix( 0, 2 ) = conic[ 3 ] * static_cast< T > ( 0.5 );
		matrix( 2, 1 ) = matrix( 1, 2 ) = conic[ 4 ] * static_cast< T > ( 0.5 );
		matrix( 2, 2 ) = conic[ 5 ];
		return Math::Matrix< T, 3, 3 >( matrix );
	}
};

/**
 * @ingroup math
 * Changes the representation of a conic from matrix to vectorial.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicFromMatrix
	: public std::unary_function< Math::Matrix< T, 3, 3 >, Math::Vector< T, 6 > >
{
public:

	/**
	 * @ingroup math
	 * Changes the representation of a conic from vectorial to matrix.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param matrix the 3x3-matrix representation of the conic
	 * @return the conic as a 6-vector
	 */
	Math::Vector< T, 6 > operator() ( const Math::Matrix< T, 3, 3 >  &matrix ) const
	{
		Math::Vector< T, 6 > conic;
		conic[ 0 ] = matrix( 0, 0 );
		conic[ 1 ] = ( matrix( 1, 0 ) + matrix( 0, 1 ) );
		conic[ 2 ] = matrix( 1, 1 );
		conic[ 3 ] = ( matrix( 2, 0 ) + matrix( 0, 2 ) );
		conic[ 4 ] = ( matrix( 2, 1 ) + matrix( 1, 2 ) );
		conic[ 5 ] = matrix( 2, 2 );
		return conic;
	}
};

/**
 * @ingroup math
 * Inverts a given conic, which usually transforms a point-conic
 * into a line-conic and vice versa due to duality.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicInverse
	: public std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 6 > >
{
public:

	/**
	 * @ingroup math
	 * Inverts a given conic, such that due to duality point-conic can be
	 * transformed into point-conics and vice versa.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the conic to be inverted
	 * @return the dual (=inverse) conic
	 */
	Math::Vector< T, 6 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T a = conic( 0 );
		const T b = conic( 1 );
		const T c = conic( 2 );
		const T d = conic( 3 );
		const T e = conic( 4 );
		const T f = conic( 5 );
		const T divisor =  static_cast< T >( 1 ) /( (a*(e*e)+c*(d*d)+(b*b)*f-a*c*f*4-b*d*e) );
		
		Math::Vector< T, 6 > iConic;
		iConic( 0 ) = -(c*f*4-e*e) * divisor;
		iConic( 1 ) =  2 * (b*f*2-d*e) * divisor;
		iConic( 2 )  = -(a*f*4-d*d) * divisor;
		iConic( 3 )  = 2 * (-(b*e-c*d*2) * divisor );
		iConic( 4 )  = 2 * (a*e*2-b*d) * divisor;
		iConic( 5 )  = -(a*c*4-b*b) * divisor;
		return Math::Vector< T, 6 >( iConic );
	};
};


/**
 * @ingroup math
 * Determines the determinant of a conic.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicDeterminant
	: public std::unary_function< Math::Vector< T, 6 >, T >
{
public:
	/**
	 * @ingroup math
	 * Calculates the determinant of a conic, based on the matrix notation.
	 *
	 * @param a 6-vector describing the conic parameters
	 * @return the determinant of the conic
	 */
	T operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T a = conic( 0 );
		const T b = conic( 1 );
		const T c = conic( 2 );
		const T d = conic( 3 );
		const T e = conic( 4 );
		const T f = conic( 5 );
		return a*c*f + (-b*b*f + b*e*d - c*d*d - a*e*e) * static_cast< T >( 0.25 );
	}
};


/**
 * @ingroup math
 * Determines the angle of a given conic. The angle expresses the
 * angular relationship between the x-axis and the major semi-axis
 * of the conic.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicAngle
	: public std::unary_function< Math::Vector< T, 6 >, T >
{
public:
	/**
	 * @ingroup math
	 * Determines the angular relationship between the
	 * x-axis and a conics' major semi-axis.
	 *
	 * The code is based on the information from the
	 * following website http://members.chello.at/gut.jutta.gerhard/kegelschnitte9.htm .
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic.
	 * @return the angle
	 */
	T operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T a = conic( 0 );
		const T b = conic( 1 );
		const T c = conic( 2 );
		const T angle = std::atan( b  / ( a - c ) ) * static_cast< T > ( 0.5 ); 
		// return 0.5 * std::atan( ( 2.0 * b ) / ( a - c ) ); 
		return ( a <= c ? angle : ( static_cast< T > ( M_PI * static_cast< T > ( 0.5 ) ) ) + angle );
	}
};

/**
 * @ingroup math
 * Determines the size(=length) of a conic's semi-axes.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicSemiAxes
	: public std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 2 > >
{
protected:
	ConicAngle< T > m_angulator;
public:

	/**
	 * @ingroup math
	 * Determines the size(=length) of a conics' semi-axes.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic.
	 * @return the size of the conics' semi-axes.
	 */
	Math::Vector< T, 2 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T theta0 = m_angulator( conic );
		const T a = conic( 0 );
		const T b = conic( 1 );// * static_cast< T > ( 0.5 );
		const T c = conic( 2 );
		const T d = conic( 3 );// * static_cast< T > ( 0.5 );
		const T e = conic( 4 );// * static_cast< T > ( 0.5 );
		const T f = conic( 5 );

		// computation shortcuts
		const T cot = std::cos(theta0);
		const T sit = std::sin(theta0);
		// double co2t = std::cos(2*theta0);
		// double si2t = std::sin(2*theta0);
		const T cot2 = cot*cot;
		const T sit2 = sit*sit;

		// Compute coefficients of the conic rotated around origin
		const T a1 = a*cot2+b*sit*cot+c*sit2;
		// double b1 = si2t*(c-a)+b*co2t; // should be equal to zero
		const T c1 = a*sit2-b*sit*cot+c*cot2;
		const T d1 = d*cot+e*sit;
		const T e1 = -d*sit+e*cot;
		const T f1 = f;

		const T num = (c1*d1*d1+a1*e1*e1-4*a1*c1*f1)/(4*a1*c1);
		const T at = std::sqrt( num/a1 );
		const T bt = std::sqrt( num/c1 );
	
		return Math::Vector< T, 2 >( at,bt );
	}

	/**
	 * @ingroup math
	 * Determines the size(=length) of a conics' semi-axes.
	 *
	 * The algorithm is based on formulas 21 and 22 on
	 * http://mathworld.wolfram.com/Ellipse.html .
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic.
	 * @return the size of the conics' semi-axes.
	 */
	/*Math::Vector< T, 2 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T a = conic( 0 );
		const T b = conic( 1 );// * static_cast< T > ( 0.5 );
		const T c = conic( 2 );
		const T d = conic( 3 );// * static_cast< T > ( 0.5 );
		const T e = conic( 4 );// * static_cast< T > ( 0.5 );
		const T f = conic( 5 );		
		const T upper = 2.0 * ( a*e*e + c*d*d + f*b*b - b*d*e - a*c*f );
		const T lower1 = b*b - a*c;
		const T lower2 = std::sqrt( std::pow( a - c, 2.0 ) + 4.0*b*b );
		const T lower3 = ( a + c ) ;
		const T aPrime = std::sqrt( upper / ( lower1 * ( lower2 - lower3 )  ) );
		const T bPrime = std::sqrt( upper / ( lower1 * ( - lower2 - lower3 )  ) );
		return Math::Vector< T, 2 >( aPrime, bPrime );
	}*/
};


/**
 * @ingroup math
 * Determines the center of a given ellipse.
 *
 * This functor should only be applied to ellipses,
 * although it is named as conic_center.
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicCenter
	: public std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 2 > >
{
public:
	/**
	 * @ingroup math
	 * Determines the center of a conic as a 2-vector.
	 *
	 * This functor should only be applied to ellipses,
	 * although it is named as conic_center.
	 * The algorithm is based on formulas 19 and 20 on
	 * http://mathworld.wolfram.com/Ellipse.html .
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic.
	 * @return the center of the conic
	 */
	Math::Vector< T, 2 > operator() ( const Math::Vector< T, 6 >& conic ) const
	{
		const T a = conic( 0 );
		const T b = conic( 1 ) * static_cast< T > ( 0.5 );
		const T c = conic( 2 );
		const T d = conic( 3 ) * static_cast< T > ( 0.5 );
		const T e = conic( 4 ) * static_cast< T > ( 0.5 );
		// const T f = conic( 5 );
	
		// (b^2)-ac
		const T divisor = static_cast< T >( 1 ) / ( b*b - a*c );
		if( divisor == 0 )
			UBITRACK_THROW( "Could not calculate the center, divisor equals zero" );
		// c*d-b*e / divisor
		const T x = ( c*d - b*e ) * divisor;
		// a*e-b*d / divisor
		const T y = ( a*e - b*d ) * divisor;
		return Math::Vector< T, 2 >( x, y );
	}
};


/**
 * @ingroup math
 * Determines the eccentricity of a given conic.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicEccentricity
	: public std::unary_function< Math::Vector< T, 6 >, T >
{
protected:
	/// helper functor to estimate the determinant
	const ConicDeterminant< T > m_determiner;
public:

	/** Standard constructor */
	ConicEccentricity( ) :
		std::unary_function< Math::Vector< T, 6 >, T >()
		, m_determiner( )
		{};
	/**
	 * @ingroup math
	 * Calculates the eccentricity of a conic.
	 *
	 * @param input conic
	 * @return the eccentricity
	 */
	T operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T a = conic( 0 );
		const T b = conic( 1 ) * static_cast< T > ( 0.5 );
		const T c = conic( 2 );
		// const T d = conic( 3 ) * static_cast< T > ( 0.5 );
		// const T e = conic( 4 ) * static_cast< T > ( 0.5 );
		// const T f = conic( 5 );
		const T upper = std::sqrt( std::pow( a-c, 2 ) + b*b );
		const T ac = (a+c);
		const T det = m_determiner( conic );
		return ( det < 0 ) ? std::sqrt( ( 2 * upper ) / (upper+ac) ) : std::sqrt( ( 2 * upper) / (upper-(ac) ) );
	}
};


/**
 * @ingroup math
 * Estimates the area of a conic.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicArea
	: public std::unary_function< Math::Vector< T, 6 >, T >
{
protected:
	const Geometry::ConicSemiAxes< T > m_axerer;
	
public:
	/** Standard constructor */
	ConicArea( ) :
		std::unary_function< Math::Vector< T, 6 >, T >()
		, m_axerer( )
		{};
	/**
	 * @ingroup math
	 * Calculates the area of a conic.
	 * 
	 *
	 * @param conic input conic
	 * @return the area of the conic
	 */
	T operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const Math::Vector< T, 2 > semi_axes = m_axerer( conic );
		return ( static_cast< T > ( M_PI ) * semi_axes( 0 ) * semi_axes( 1 ) );
	};
};


/**
 * @ingroup math
 * Signs if conic is a circle.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct IsConicCircle
	: public std::unary_function< Math::Vector< T, 6 >, bool >
{
protected:
	const T m_error;
public:
	/**
	 * @ingroup math
	 * Constructor call to set an epsilon for the decision.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param epsilon the minimal error to be tolerated for b being zero
	 */
	IsConicCircle( const T error = static_cast< T > ( 1e-3 ) )
		: std::unary_function< Math::Vector< T, 6 >, bool >( )
		, m_error( error )
	{ }

	/**
	 * @ingroup math
	 * Signs if conic is a circle.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic
	 * @return flag if conic is circle or not 
	 */
	bool operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		// b ~ 0
		if( std::fabs( conic( 1 ) ) > m_error )
			return false;
		// |a-c| ~ 0
		if( std::fabs( conic( 0 ) - conic( 2 ) ) > m_error  )
			return false;
		
		//bb-4ac < 0
		return( ( conic( 1 ) *  conic( 1 ) - 4 *  conic( 0 ) *  conic( 2 ) ) < 0 );
	}
};


/**
 * @ingroup math
 * Signs if conic is degenerate.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct IsConicDegenerate
	: public std::unary_function< Math::Vector< T, 6 >, bool >
{
protected:

	const T m_epsilon;
	const Geometry::ConicDeterminant< T > m_determiner;
	
public:
	/**
	 * @ingroup math
	 * Constructor call to set an epsilon for the decision.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param epsilon the minimal error to be tolerated for the determinant
	 */
	IsConicDegenerate( const T epsilon = static_cast< T > ( 1e-3 ) )
		: std::unary_function< Math::Vector< T, 6 >, bool >()
		, m_epsilon( epsilon )
		, m_determiner( )
		{}
	
	/**
	 * @ingroup math
	 * Calculates if a conic is degenerated, by calculating the determinant.
	 *
	 * @param conic the input conic
	 * @return signs if conic is degenerate or not
	 */
	bool operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		return std::fabs( m_determiner( conic ) ) < m_epsilon ;
	}
};

/**
 * @ingroup math
 * Signs if conic is an ellipse.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct IsConicEllipse
	: public std::unary_function< Math::Vector< T, 6 >, bool >
{
public:
	/**
	 * @ingroup math
	 * Signs if conic is a ellipse.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic
	 * @return flag if conic is ellipse or not 
	 */
	bool operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		//ellipses must satisfy : bb-4ac < 0 
		return ( conic( 1 ) * conic( 1 ) - 4 * conic( 0 ) * conic( 2 ) ) < 0 ? true : false;
	}
};

/**
 * @ingroup math
 * Signs if conic is an hyperbola.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct IsConicHyperbola
	: public std::unary_function< Math::Vector< T, 6 >, bool >
{
public:
	/**
	 * @ingroup math
	 * Signs if conic is an hyperbola.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic
	 * @return flag if conic is hyperbola or not 
	 */
	bool operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		//ellipses must satisfy : bb-4ac > 0 
		return ( conic( 1 ) * conic( 1 ) - 4 * conic( 0 ) * conic( 2 ) ) > 0 ? true : false;
	}
};


/**
 * @ingroup math
 * Signs if conic is a parabola.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct IsConicParabola
	: public std::unary_function< Math::Vector< T, 6 >, bool >
{
protected:
	const T m_epsilon;
	
public:
	/**
	 * @ingroup math
	 * Constructor call to set an epsilon for the decision.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param epsilon tolerance value for accepting a conic as parabola
	 */
	IsConicParabola( const T epsilon = 1e-2 )
		: std::unary_function< Math::Vector< T, 6 >, bool >()
		, m_epsilon( epsilon )
		{}
	
	/**
	 * @ingroup math
	 * Signs if conic is a parabola.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic
	 * @return flag if conic is parabola or not 
	 */
	bool operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		// parabola satisfy : bb-4ac == 0 
		// epsilon is applied since it might never be exactly zero
		return std::fabs( conic( 1 ) * conic( 1 ) - 4 * conic( 0 ) * conic( 2 ) ) < m_epsilon;
	}
};


/**
 * @ingroup math
 * Scales the semi-axes of a given conic.
 * 
 * Use this function carefully since it can alter
 * the conics' position in space.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ScaleConicUnsafe
	: public std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 6 > >
{
protected:
	/// scale for major semi-axes
	const T m_scaleA1;
	/// quadratic scale for major semi-axes
	const T m_scaleA2;
	/// scale for minor semi-axes
	const T m_scaleB1;
	/// quadratic scale for minor semi-axes
	const T m_scaleB2;

public:
	/**
	 * @ingroup math
	 * Constructor call to set the scaling parameters for all conics.
	 *
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param scaleA the scale for the major semi-axis
	 * @param scaleB the scale for the minor semi-axis
	 */
	ScaleConicUnsafe( const T scaleA, const T scaleB )
		: std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 6 > >( )
		, m_scaleA1(  static_cast< T >( 1 ) / scaleA )
		, m_scaleA2(  m_scaleA1 * m_scaleA1 )
		, m_scaleB1(  static_cast< T >( 1 ) / scaleB )
		, m_scaleB2(  m_scaleB1 * m_scaleB1 )
	{
		// just for illustration what happens here
		// [  C1_1*sa^2, C1_2*sa*sb, C1_3*sa]
		// [ C2_1*sa*sb,  C2_2*sb^2, C2_3*sb]
		// [    C3_1*sa,    C3_2*sb,    C3_3]
	};
	
	/**
	 * @ingroup math
	 * Constructor call to set the scaling parameters fo all conics.
	 *
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param scale the scale for the both the semi-axes
	 */
	ScaleConicUnsafe( const T scale )
		: std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 6 > >( )
		, m_scaleA1( static_cast< T >( 1 ) / scale )
		, m_scaleA2(  m_scaleA1 * m_scaleA1 )
		, m_scaleB1(  m_scaleA1 )
		, m_scaleB2(  m_scaleA2 )
	{
		// just for illustration what happens here		
		// [ a/s1^2, b/s1^2, d/s1]
		// [ b/s1^2, c/s1^2, e/s1]
		// [   d/s1,   e/s1,    f]
	};

	/**
	 * @ingroup math
	 * Scales the given conic by the scales handed over in the
	 * constructor.
	 *
	 * Attention: This scaling can alter the conics' position
	 * in space.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic the input conic.
	 * @return the scaled conic. 
	 */
	Math::Vector< T, 6 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		Math::Vector< T, 6 > vec;
		const T a = conic[ 0 ];
		const T b = conic[ 1 ];
		const T c = conic[ 2 ];
		const T d = conic[ 3 ];
		const T e = conic[ 4 ];
		const T f = conic[ 5 ];
		vec[ 0 ] = a * m_scaleA2;
		vec[ 1 ] = b * m_scaleA1*m_scaleB1;
		vec[ 2 ] = c * m_scaleB2;
		vec[ 3 ] = d * m_scaleA1;
		vec[ 4 ] = e * m_scaleB1;
		vec[ 5 ] = f;
		return Math::Vector< T, 6 >( vec );
	}
};


/**
 * @ingroup math
 * Translates a conic relatively in space by a given 2-vector.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct TranslateConic
	: public std::binary_function< Math::Vector< T, 2 >, Math::Vector< T, 6 >, Math::Vector< T, 6 > >
{
public:
	/**
	 * @ingroup math
	 * Translates a conic relatively in space by the given 2-vector.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param translation expresses the relative positioning
	 * @param conic the input conic.
	 * @return the translated conic. 
	 */
	Math::Vector< T, 6 > operator() ( const Math::Vector< T, 2 > &translation, const Math::Vector< T, 6 > &conic ) const
	{
		
		const T a = conic( 0 );
		const T b = conic( 1 ) * static_cast< T > ( 0.5 );
		const T c = conic( 2 );
		const T d = conic( 3 ) * static_cast< T > ( 0.5 );
		const T e = conic( 4 ) * static_cast< T > ( 0.5 );
		const T f = conic( 5 );
		const T tx ( translation( 0 ) );
		const T ty ( translation( 1 ) );
		const T atx = a * tx;
		const T aty = a * ty;
		const T btx = b * tx;
		const T bty = b * ty;
		
		// this is the matrix expression for the translation
		// T[0][0] = a;
		// T[0][1] = b;
		// T[0][2] = d-a*a31-a32*b;
		// T[1][0] = b;
		// T[1][1] = c;
		// T[1][2] = e-a31*b-a32*c;
		// T[2][0] = d-a*a31-a32*b;
		// T[2][1] = e-a31*b-a32*c;
		// T[2][2] = f-a31*d-a32*e+a31*(-d+a*a31+a32*b)+a32*(-e+a31*b+a32*c);
		
		Math::Vector< T, 6 > vec;
		vec( 0 ) = conic( 0 );
		vec( 1 ) = conic( 1 );
		vec( 2 ) = conic( 2 );

		vec( 3 ) = d - atx - bty;
		vec( 4 ) = e - btx - c*ty;
		vec( 5 ) = f - d*tx - e*ty + tx*(-d + atx + bty) + ty*(-e + btx + ty*c);
		
		vec( 3 ) *= 2;
		vec( 4 ) *= 2;
		return Math::Vector< T, 6 >( vec );
	}
};

/**
 * @ingroup math
 * Determines a  p^T *C* p product from a given conic C and a
 * pixel p.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicPixel
	: public std::binary_function< Math::Vector< T, 6 >, Math::Vector< T, 2 >, T >
{
public:
	/**
	 * @ingroup math
	 * Determines a  vector matrix vector product (=v^tCv)
	 * from a given conic C and a pixel p.
	 *
	 * This function can be used to determine whether a pixel is
	 * on the outline of a conic or inside or outside a conic.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param conic input conic
	 * @param pixel input pixel
	 * @return resulting product
	 */
	T operator() ( const Math::Vector< T, 6 > &conic, const Math::Vector< T, 2 > &pixel  ) const
	{
		const T x = pixel[ 0 ];
		const T y = pixel[ 1 ];
		
		const T a = conic[ 0 ];
		const T b = conic[ 1 ] * static_cast< T > ( 0.5 );
		const T c = conic[ 2 ];
		const T d = conic[ 3 ] * static_cast< T > ( 0.5 );
		const T e = conic[ 4 ] * static_cast< T > ( 0.5 );
		const T f = conic[ 5 ];
		
		return ( f + e*y + d*x + x*( d + a*x + b*y) + y *( e + b*x + c*y ) );
	}
};

/**
 * @ingroup math
 * Reflects a conic at a given y-axis.
 * 
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct FlipConicHorizontal
	: public std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 6 > >
{
protected:
	/** signs where to reflect the conic at the y-axis */
	const T m_y;
	
public:

	/**
	 * @ingroup math
	 * Constructor call to set the y-axis at which a conic will be
	 * reflected.
	 *
	 * This function can be used to flip a conic estimated in an image
	 * with the origin flag set to TOPLEFT to a conic in an image 
	 * with origin at BOTTOMLEFT.
	 *
	 * @tparam T type of conic ( e.g. \c double or \c float )
	 * @param height the value of the y-axis, usually image height (will be subtracted by one)
	 */
	FlipConicHorizontal( const T height )
		: std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 6 > >( )
		, m_y( height - 1  )
	{ }
	
	Math::Vector< T, 6 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T a = conic[ 0 ];
		const T b = conic[ 1 ] * static_cast< T > ( 0.5 );
		const T c = conic[ 2 ];
		const T d = conic[ 3 ] * static_cast< T > ( 0.5 );
		const T e = conic[ 4 ] * static_cast< T > ( 0.5 );
		const T f = conic[ 5 ];
		
		// the following matrix expresses the conic flipping.
		// [       a,        -b,               d + b*y]
		// [      -b,         c,             - e - c*y]
		// [ d + b*y, - e - c*y, f + e*y + y*(e + c*y)]
		Math::Vector< T, 6 > vec;
		vec[ 0 ] = a;
		vec[ 1 ] = -b;
		vec[ 2 ] = c;
		vec[ 3 ] = d + b * m_y;
		vec[ 4 ] = -( e + c * m_y );
		vec[ 5 ] = f +e*m_y + m_y*( e+c*m_y);
		vec[ 1 ] *= 2;
		vec[ 3 ] *= 2;
		vec[ 4 ] *= 2;
		return Math::Vector< T, 6 >( vec );
	}
};


/**
 * @ingroup math
 * Determines the upper and lower limits (y-axis) of a conic.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicUpperLowerLimit
	: public std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 2 > >
{
public:
	/**
	 * @ingroup math
	 * Determines the upper and lower limits (y-axis) of a conic.
	 *
	 * Values are returned in the order lower than upper limit
	 * in a positive direction of the y-axis. This might result in
	 * different order for images with their origin at the upper left.
	 *
	 * @param conic input conic
	 * @return the upper and lower limits of the conic
	 */
	Math::Vector< T, 2 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{

		const T a = conic( 0 );
		const T b = conic( 1 ) * static_cast< T > ( 0.5 );
		const T c = conic( 2 );
		const T d = conic( 3 ) * static_cast< T > ( 0.5 );
		const T e = conic( 4 ) * static_cast< T > ( 0.5 );
		const T f = conic( 5 );

		T lower = 4 * ( b*b - a*c );
		const T upper1 = 8 * ( b*d - a*e );
		const T temp = 4 * ( d*d - a*f );
		const T upper2 = std::sqrt( upper1 * upper1 - ( 4 * lower * temp ) );
		lower *= 2;
		const T y1 = -( (upper1/lower) + (upper2/lower) );
		const T y2 = -( (upper1/lower) - (upper2/lower) );
	
		if ( y1 < y2 )
			return Math::Vector< T, 2 >( y1, y2 );
		else
			return Math::Vector< T, 2 >( y2, y1 );
	}
};


/**
 * @ingroup math
 * Determines the left and right limit (x-axis) of a conic.
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicLeftRightLimit
	: public std::unary_function< Math::Vector< T, 6 >, Math::Vector< T, 2 > >
{
public:
	/**
	 * @ingroup math
	 * Determines the left and right limit (x-axis) of a conic.
	 *
	 * Values are returned in the order left than right in respect
	 * to the x-axis.
	 *
	 * @param conic input conic
	 * @return the left and right limits of the conic
	 */
	Math::Vector< T, 2 > operator() ( const Math::Vector< T, 6 > &conic ) const
	{
		const T a = conic( 0 );
		const T b = conic( 1 ) * static_cast< T > ( 0.5 );
		const T c = conic( 2 );
		const T d = conic( 3 ) * static_cast< T > ( 0.5 );
		const T e = conic( 4 ) * static_cast< T > ( 0.5 );
		const T f = conic( 5 );
	
		const T x1 = (b*e-c*d+c*std::sqrt((a*(e*e)+c*(d*d)+(b*b)*f-a*c*f-b*d*e*2)/c))/(a*c-b*b);
		const T x2 = -(-b*e+c*d+c*std::sqrt((a*(e*e)+c*(d*d)+(b*b)*f-a*c*f-b*d*e*2)/c))/(a*c-b*b);

		return ( x1 < x2 ) ? Math::Vector< T, 2 >( x1, x2 ) : Math::Vector< T, 2 >( x2, x1 );
    }
};


/**
 * @ingroup math
 * Determines left/right x-values of a conic from the 
 * intersection of the conic and a line parallel to the y-axis (0 slope).
 *
 * This functor class can be applied to STL-containers
 * of conics via a STL-algorithms.
 */
template< typename T >
struct ConicHorizontalIntersection
	: public std::binary_function< Math::Vector< T, 6 >, T, Math::Vector< T, 2 > >
{
public:
	/**
	 * @ingroup math
	 * Creates a quadric-matrix from a given position and their semi-axes
	 *
	 * @param conic the input conic
	 * @param y the height of the line
	 * @return the left/right x-values of the intersection
	 */
	Math::Vector< T, 2 > operator() ( const Math::Vector< T, 6 > &conic, const T y ) const
	{
		const T b = conic( 1 )*0.5 * y + conic( 3 ) * static_cast< T > ( 0.5 ); 
		const T c = ( conic( 2 ) * y + conic( 4 ) ) * y + conic( 5 );
		const T d = std::sqrt( b * b - conic( 0 ) * c );
		const T x1 = ( (-b + d) / conic( 0 ) );
		const T x2 = ( (-b - d) / conic( 0 ) );
		return ( x1 < x2 ) ? Math::Vector< T, 2 >( x1, x2 ) : Math::Vector< T, 2 >( x2, x1 );
	}
};


} } } // namespace Ubitrack::Math::Geometry

#endif  // __H__CONIC_FUNCTORS__

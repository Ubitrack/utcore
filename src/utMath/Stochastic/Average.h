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
 * @ingroup dataflow_components
 * @file
 * This class calculates the average of a List measurement.
 *
 * @author Florian Echtler <echtler@in.tum.de>
 * @author Christian Waechter <christian.waechter@in.tum.de> (modified)
 */

#ifndef __UBITRACK_MATH_STOCHASTICAVERAGE_INCLUDED__
#define __UBITRACK_MATH_STOCHASTICAVERAGE_INCLUDED__

// Ubitrack
#include "../Scalar.h"
#include "../Pose.h"
#include "../ErrorPose.h"
#include "../Vector.h"
#include "../ErrorVector.h"

#include "../Blas2.h" // outer_product
#include "../Util/TypeToVector.h" // castToVector

namespace Ubitrack { namespace Math { namespace Stochastic {

/**
 * @brief Data structure to support online averaging of measurements.
 *
 * This struct can be used to calculate the mean value of any known relevant measurement type
 * There are several specializations of this struct for special purposes like
 * calculating a covariance as well or providing a mean for built-in types.
 * 
 * The struct can be applied to sequence container in the simplest way using \c std::for_each .
 * Using \c std::for_each the Average can be applied on various measurement types in various kinds
 * of sequence containers (e.g. \c std::vector, \c std::list or \c std::set, etc ).
 *
 * The following calculations are supported at the moment:
 * - \b Scalar< T > -> \b Scalar< T >
 * - \b Quaternion -> \b Quaternion
 * - \b Vector< T, N > -> \b Vector< T, N >
 * - \b Vector< T, N > -> \b ErrorVector< T, N >
 * - \b Pose -> \b Pose
 * - \b Pose -> \b ErrorPose
 * - built-in (e.g. \c double ) -> built-in (e.g. \c float )
 *
 * Example use case illustrating the struct :\n
 @code
 #include <algorithm> // std::for_each
 
 std::vector< Vector3d > points3d; // <- should be filled with values
 Average< Vector3d > Averager;
 Averager = std::for_each( points3d.begin(), points3d.end(), Averager );
 Averager.getAverage();
 
 // or
 Average< ErrorVector< double, 3 > > AveragerCovariance;
 AveragerCovariance = std::for_each( points3d.begin(), points3d.end(), AveragerCovariance );
 AveragerCovariance.getAverage();
 
 @endcode
 * 
 */
template< typename ResultType, std::size_t N = Util::TypeToVector< ResultType >::size >
struct Average
{
	/** defines the type that should be used as output type for the measurement. */
	typedef ResultType value_type;
	
	/** defines the precision of the underlying built-in type (e.g. \c double or \c float )*/
	typedef typename Util::TypeToVector< value_type >::precision_type precision_type;
	
	/** defines the vector type used to summarize all values*/
	typedef typename Math::Vector< precision_type, N > mean_type;
	
	/** Keeps the count of the amount of elements already pushed into the Average */
	std::size_t m_counter;
	
	/** The vector containing the sums of all received measurement so far. */
	mean_type m_mean;
	
	/** Standard constructor */
	Average()
		: m_counter( 0 )
		, m_mean ( mean_type::zeros() )
	{}
	
	/** unary bracket operator accepting measurements that are added to the mean value. */
	void operator() ( const value_type& value )
	{
		++m_counter;
		mean_type tmp;
		Util::castToVector( value, tmp );
		m_mean += tmp;
	}
	
	/** function that calculates the mean value and returns it.*/
	value_type getAverage() const
	{
		return m_mean / m_counter;
	}
};

/// @internal specialization for built-in types and/or Math::Scalar
template< typename ResultType >
struct Average< ResultType, 1 >
{
	typedef ResultType value_type;
	typedef typename Util::TypeToVector< value_type >::precision_type precision_type;
	typedef precision_type mean_type;
	
	std::size_t m_counter;
	mean_type m_mean;
	
	Average()
		: m_counter( 0 )
		, m_mean ( 0 )
	{}
	
	template< typename T >
	void operator() ( const T& value )
	{
		++m_counter;
		m_mean += static_cast< value_type > ( value );
	}
	
	value_type getAverage() const
	{
		return static_cast< value_type > ( m_mean / m_counter );
	}
};

/// overloaded unary bracket operator for Quaternion measurements, that brings all rotation measurements in the same hemisphere
template<>
void Average< Math::Quaternion >::operator() ( const value_type& value )
{
	++m_counter;
	mean_type tmp;
	Util::castToVector( value, tmp );
	m_mean += (( value.w() >= 0 ) ? tmp : tmp * (-1) );
};

/// overloaded getAverage function for Quaternion measurements, when new struct is available this should not be necessary anymore.
template<>
Math::Quaternion Average< Math::Quaternion >::getAverage() const
{
	const mean_type mean = m_mean / m_counter;
	return Quaternion( mean[ 0 ], mean[ 1 ], mean[ 2 ], mean[ 3 ] ).normalize();
};

/// overloaded unary bracket operator for pose measurements, that brings all rotation measurements in the same hemisphere
template<>
void Average< Math::Pose >::operator() ( const value_type& value )
{
	++m_counter;
	mean_type tmp;
	Util::castToVector( value, tmp );
	// The order is tx, ty, tz, qx, qy, qz, qw.
	// check now for real part (==qw) being positive and correct if necessary:
	if( tmp[ 6 ] < 0 )
	{
		tmp[ 3 ] *= (-1);
		tmp[ 4 ] *= (-1);
		tmp[ 5 ] *= (-1);
		tmp[ 6 ] *= (-1);
	}
	m_mean += tmp;
};

/// overloaded getAverage function for Pose measurements, when new struct is available this should not be necessary anymore as well.
template<>
Math::Pose Average< Math::Pose >::getAverage() const
{
	
	const mean_type mean = m_mean / m_counter;
	const Math::Pose pose = Math::Pose::fromVector( mean );
	return Math::Pose( Math::Quaternion( pose.rotation() ).normalize(), pose.translation() );
};

/// @internal specialization of average struct for an ErrorVector measurement as mean+covariance.
template< typename T, std::size_t N >
struct Average< Math::ErrorVector< T, N >, N >
	: public Average< Math::Vector< T, N > >
{
private:
	typedef Average< Math::Vector< T, N > > super;
	
public:
	typedef Math::ErrorVector< T, N > value_type;
	typedef Math::Matrix< T, N, N > varianz_type;
	
	varianz_type m_covariance;
	
	Average()
		: super()
		, m_covariance ( varianz_type::zeros() )
	{}
	
	void operator() ( const typename super::value_type& value )
	{
		super::operator( )( value );
		m_covariance += Math::outer_product( value, value );
	}

	value_type getAverage() const
	{
		const typename super::value_type mean = super::getAverage();
		const varianz_type covarianz = ( m_covariance / super::m_counter );
		return value_type( mean, covarianz - Math::outer_product( mean, mean ) );
	};
};

/// @internal specialization of average struct for an ErrorPose measurement as mean+covariance.
template<>
struct Average< Math::ErrorPose >
	: public Average< Math::Pose >
{
private:
	typedef Average< Math::Pose > super;
	
public:
	typedef Math::ErrorPose value_type;
	typedef value_type::value_type precision_type;
	static const Util::TypeToVector< super::value_type >::size_type size = Util::TypeToVector< super::value_type >::size;
	typedef Math::Matrix< precision_type, size, size > varianz_type;
	
	varianz_type m_covariance;
	
	Average()
		: super()
		, m_covariance ( varianz_type::zeros() )
	{}
	
	void operator() ( const super::value_type& value )
	{
		Math::Vector< precision_type, 7 > tmp;
		Util::castToVector( value, tmp );
		
		// The order is tx, ty, tz, qx, qy, qz, qw.
		// check now for real part (==qw) being positive and correct if necessary:
		if( tmp[ 6 ] < 0 )
		{
			tmp[ 3 ] *= (-1);
			tmp[ 4 ] *= (-1);
			tmp[ 5 ] *= (-1);
			tmp[ 6 ] *= (-1);
		}
		++m_counter;
		m_mean += tmp;
		m_covariance += Math::outer_product( tmp, tmp );
	}

	value_type getAverage() const
	{
	
		/*
		 * Use inverted mean value to transform the additive 7x7
		 * covariance to the 6x6 multiplicative format The conversion is
		 * conducted according to the following formulas:
		 * 
		 * q_m = q_0 * ( q_id + q_e )
		 * 
		 * where q_id is the identity quaternion and q_e is a quaternion
		 * with expectation ((0,0,0),0) and a covariance covering only the
		 * imaginary part. Together ( q_id + q_e ) represent a small
		 * quaternion ((e_rx, e_ry, e_rz), 1). If mean and covariance of
		 * the quaternion are estimated according to the usual formulas,
		 * however, one gets the following instead:
		 * 
		 * q_m = q_0 + q'_e
		 * 
		 * Together with the first formula, this yields
		 * 
		 * q_0 * ( q_id + q_e ) = q_0 + q'_e
		 * ( q_id + q_e )       = ~q_0 * q_0 + ~q_0 * q'_e
		 * q_e                  = q_id + ~q_0 * q'_e - q_id
		 * q_e                  = ~q_0 * q'_e
		 *
		 * Thus, one has to rotate the distribution by ~q_0. The variance
		 * of the real part can then be discarded, it should be ~0.
		 */
		const super::value_type meanPose = super::getAverage();
		
		super::mean_type meanVec;
		Util::castToVector( meanPose, meanVec );
		
		super::mean_type invMeanVec;
		Util::castToVector( (~meanPose), invMeanVec );
		
		
		const varianz_type covarianz = ( m_covariance / m_counter );
		Math::ErrorVector< precision_type, 7 > ev( invMeanVec, covarianz - Math::outer_product( meanVec, meanVec ) );
		const Math::ErrorPose invEp = Math::ErrorPose::fromAdditiveErrorVector( ev );
		// We created the error pose from the inverted mean value above, to obtain the transformed 6x6 covariance
		// Now, we recreate the error pose with the computed mean value.	
		return value_type( meanPose, invEp.covariance() );
	};
	
};

} } } // namespace Ubitrack::Math::Stochastic
#endif //__UBITRACK_MATH_STOCHASTICAVERAGE_INCLUDED__
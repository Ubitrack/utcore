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

// Ubitrack
#include <utUtil/Exception.h>

// Boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/constants/constants.hpp>

#include <utMath/Pose.h>
#include <utMath/ErrorPose.h>
#include <utMath/Vector.h>
#include <utMath/ErrorVector.h>

namespace Ubitrack { namespace Tracking {

namespace ublas = boost::numeric::ublas;

template< class EventType, class ResultType >
class Average
{
public:
	/**
	 * Standard component constructor.
	 *
	 * @param nm Unique name of the component.
	 * @param cfg ComponentConfiguration containing all configuration.
	 */
	Average(  )		
	{
	}

	typedef typename std::vector< EventType > EventList;
	
	ResultType mean( const EventList &eList )
	{
		ResultType tmp;
		size_t size = eList.size();
		
		for( size_t i = 0; i < size; i++ )
			tmp = tmp + ( eList[i] / size);
		return tmp;
	};

protected:

	

	
	ResultType incrementalEstimate( EventType& perturbed );
		
	ublas::vector< double > meanv;
	ublas::matrix< double > outProd;

	int m_counter;
};

template<>
Math::Vector< 3 > Average< Math::Vector< 3 >, Math::Vector< 3 > >::mean( const std::vector< Math::Vector< 3 > > &eList )
{
	size_t size = eList.size();
	Math::Vector< 3 > m_mean ( 0.0, 0.0, 0.0 );

	for ( size_t i = 0; i < size; i++ )
	{
		m_mean += eList[i];
	}
 	return m_mean / static_cast< double >( size );
}

template<>
Math::ErrorVector< 3 > Average< Math::Vector< 3 >, Math::ErrorVector< 3 > >::mean( const std::vector< Math::Vector< 3 > > &eList )
{
	size_t size = eList.size();
	ublas::vector< double > m_mean = ublas::zero_vector< double >( 3 );
	ublas::matrix< double > m_outProd = ublas::zero_matrix< double >( 3, 3 );

	for (size_t i = 0; i < size; i++)
	{
		m_mean = m_mean + ( eList[i] / size);
		m_outProd = m_outProd + ublas::outer_prod( eList[i], eList[i] );
	}
	
	Math::ErrorVector< 3 > ev ( m_mean, m_outProd / ( static_cast< double > ( size ) ) - ublas::outer_prod ( m_mean, m_mean ) );

 	return ev;
}

template<>
Math::Quaternion Average< Math::Quaternion, Math::Quaternion >::mean( const std::vector< Math::Quaternion > &eList )
{
	size_t size = eList.size();
	Math::Quaternion q_mean ( 0.0, 0.0, 0.0, 0.0 );

	Math::Vector< 3 > axisAngle( 0.0, 0.0, 0.0 );
	
	// angle = 2 * acos( q(:,4) );
// s = sqrt( 1 - q(:,4) .* q(:,4) );% assuming quaternion normalised then w is less than 1, so term always positive.

// aa = [ angle, q(:,1:3) ];
// if( s > eps ) % test to avoid divide by zero, s is always positive due to sqrt

    // aa(:, 2) = aa(:, 2) ./ s;
    // aa(:, 3) = aa(:, 3) ./ s;
    // aa(:, 4) = aa(:, 4) ./ s;
// end
	int num( 0 );
	
	for( size_t i = 0; i < size; i++ )
	{
		Math::Quaternion q_t = eList[i];
		
		if( q_t.w() < 0.0 )
			q_t *= -1.0;
			
		double angle = 2.0 * std::acos( q_t.w() );
		double scale = std::sqrt( 1.0 - q_t.w() * q_t.w() );
		if( scale > 0.000001 )
		{
			Math::Vector< 3 > axis( q_t.x(), q_t.y(), q_t.z() );
			axis /= scale;
			axis *= angle;
			axisAngle += axis;
			++num;
		}
		
		// q_mean += q_t;
	}
	if( num > 0 )
		axisAngle /= num;
		
	double norm = ublas::norm_2( axisAngle );
	double s = std::sin( norm / 2.0 );
	axisAngle *= ( s / norm ) ;
	Math::Quaternion q_final( axisAngle( 0 ), axisAngle( 1 ), axisAngle( 2 ), std::cos( norm/ 2.0 ) );
	// q_mean /= size;
	// q_mean.normalize();
	
	return q_final;
}

template<>
Math::Pose Average< Math::Pose, Math::Pose >::mean( const std::vector< Math::Pose > &eList )
{
	size_t size = eList.size();
	Math::Vector< 3 > p_mean ( 0.0, 0.0, 0.0 );
	Math::Quaternion q_mean ( 0.0, 0.0, 0.0, 0.0 );

	for( size_t i = 0; i < size; i++ )
	{
		Math::Quaternion q_t = eList[i].rotation();
		if( q_t.w() < 0.0 )
			q_t *= -1.0;
			
		q_mean += q_t;
		p_mean += eList[i].translation();
	}
	p_mean /= (double)size;
	q_mean /= (double)size;
	q_mean.normalize();
	
	return Math::Pose( q_mean, p_mean );
}


Math::ErrorPose incEstimate(  Math::Pose poseNew,  ublas::vector< double >& meanv,  ublas::matrix< double >& outProd, int m_counter)
{
	ublas::vector_range< ublas::vector<double> > posMean( meanv, ublas::range( 0, 3 ) );
	ublas::vector_range< ublas::vector<double> > rotMean( meanv, ublas::range( 3, 7 ) );

	//LOG4CPP_TRACE ( logger, "Update pose event: " << poseNew );

	// The order is tx, ty, tz, qx, qy, qz, qw.
	ublas::vector< double > poseNewVec( 7 );
	poseNew.toVector( poseNewVec );
	ublas::vector_range< ublas::vector<double> > posNew( poseNewVec, ublas::range( 0, 3 ) );
	ublas::vector_range< ublas::vector<double> > rotNew( poseNewVec, ublas::range( 3, 7 ) );

	// Take care of quaternion ambiguity
 	if ( ublas::inner_prod( rotNew, rotMean ) < 0 )
 		rotNew *= -1;

	// Update running mean value
	meanv = ( ( ((double)m_counter - 1) / (double)m_counter ) * meanv ) + ( ( 1 / (double)m_counter ) * poseNewVec );

	// Running outer product of pose random variable (not yet normalized by number of measurements)
	outProd = outProd + ublas::outer_prod( poseNewVec, poseNewVec );

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

	Math::Vector<7> invMean;
	(~(Math::Pose::fromVector( meanv ) ) ).toVector( invMean );
	Math::ErrorVector< 7 > ev ( invMean, outProd / ( (double)m_counter ) - ublas::outer_prod ( meanv, meanv ) );
	Math::ErrorPose invEp = Math::ErrorPose::fromAdditiveErrorVector( ev );
	
	// We created the error pose from the inverted mean value above, to obtain the transformed 6x6 covariance
	// Now, we recreate the error pose with the computed mean value.
	Math::ErrorPose ep( Math::Pose::fromVector( meanv ), invEp.covariance() );

	//LOG4CPP_TRACE( logger, "Running (empirical) mean / covariance: " << std::endl << ep );

	// For debug purposes, compute positional and angular error...
	Math::Matrix< 6, 6 > covar = ep.covariance();
	double posRms = sqrt ( covar (0,0) + covar (1,1) + covar (2,2) );
	//LOG4CPP_INFO( logger, "RMS positional error [mm]: " << posRms );
	ublas::vector< double > axis (3);
	axis (0) = sqrt ( covar (3,3) );
	axis (1) = sqrt ( covar (4,4) );
	axis (2) = sqrt ( covar (5,5) );
	double norm = norm_2 (axis);
	double phi = asin ( norm ) * 2;
	phi = phi * 180 / boost::math::constants::pi<double>();
	//LOG4CPP_INFO( logger, "Standard deviation of rotational error [deg]: " << phi );
	
	return ep;
}


template<>
Math::ErrorPose Average< Math::Pose, Math::ErrorPose >::mean( const std::vector< Math::Pose > &eList )
{
	size_t size = eList.size();
	
	meanv = ublas::zero_vector< double >( 7 );
	outProd = ublas::zero_matrix< double >( 7, 7 );

	m_counter = 1;
	Math::ErrorPose estimate;
	for( size_t i = 0; i < size; i++ ) {
		
		estimate = incEstimate ( (eList.at(i)), meanv, outProd, m_counter );
		m_counter++;
	}
	return estimate;
	/*
	Math::Vector< 3 > p_mean ( 0.0, 0.0, 0.0 );
	Math::Quaternion q_mean ( 0.0, 0.0, 0.0, 0.0 );
	
	ublas::vector< double > m_mean = ublas::zero_vector< double >( 7 );
	ublas::matrix< double > m_outProd = ublas::zero_matrix< double >( 7, 7 );
	
	for( unsigned i = 0; i < size; i++ )
	{
		Math::Quaternion q_t = eList[i].rotation();
		if( q_t.w() < 0.0 )
			q_t *= -1.0;
			
		q_mean += q_t;
		p_mean += eList[i].translation();
		Math::Pose pose( q_mean, p_mean );
		
		ublas::vector< double > poseTmpVec( 7 );
		pose.toVector( poseTmpVec );
		m_outProd = m_outProd + ublas::outer_prod( poseTmpVec, poseTmpVec );
	}
	
	p_mean /= size;
	q_mean /= size;
	q_mean.normalize();
	
	Math::Pose pose_mean( q_mean, p_mean );
	
	pose_mean.toVector( m_mean );

	Math::Vector< 7 > invMean;
	(~(Math::Pose::fromVector( m_mean ) ) ).toVector( invMean );
	Math::ErrorVector< 7 > ev ( invMean, m_outProd / ( static_cast< double > ( size ) ) - ublas::outer_prod ( m_mean, m_mean ) );
	Math::ErrorPose invEp = Math::ErrorPose::fromAdditiveErrorVector( ev );
	
	Math::ErrorPose ep( Math::Pose::fromVector( m_mean ), invEp.covariance() );
	*/
	//return ep;
}


} } // namespace Ubitrack::Components


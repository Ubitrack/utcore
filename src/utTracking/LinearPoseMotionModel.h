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
 * @ingroup tracking
 * @file
 * Motion model for kalman filtering of poses
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */
 
#ifndef __UBITRACK_TRACKING_LINEARPOSEMOTIONMODEL_H_INCLUDED__
#define __UBITRACK_TRACKING_LINEARPOSEMOTIONMODEL_H_INCLUDED__


namespace Ubitrack { namespace Tracking {

/**
 * Motion model for pose with linear n-th order motions for position and orientation.
 */
class LinearPoseMotionModel
{
public:
	/**
	 * Initializes the motion model.
	 * posOrder and oriOrder set the number of derivatives for position and orientation.
	 * @param posOrder number of position derivatives. -1 means no position.
	 * @param oriOrder number or orientation derivatives. -1 means no orientation.
	 */
	LinearPoseMotionModel( int posOrder, int oriOrder )
		: m_posOrder( posOrder )
		, m_oriOrder( oriOrder )
		, m_processNoise( 3 + 3 * m_posOrder + ( m_oriOrder >= 0 ? 4 + 3 * m_oriOrder : 0 ) )
	{
	}
	
	/** returns the required size of the state vector */
	std::size_t stateSize() const
	{ return m_processNoise.size(); }
	
	/** returns the number of position derivatives */
	int posOrder() const
	{ return m_posOrder; }
	
	/** returns the number of orientation derivatives */
	int oriOrder() const
	{ return m_oriOrder; }
	
	/** 
	 * sets the process noise for position.
	 * @param order number of derivative, 0 is the absolute position
	 * @param value standard deviation per second in m/s, m/s^2, m/s^3, ..
	 */
	void setPosPN( const std::size_t order, const double value )
	{ 
		for ( std::size_t i( 3 * order ); i < (3 + 3 * order); i++ )
			m_processNoise( i ) = value * value;
	}

	/** 
	 * sets the process noise for orientation.
	 * @param order number of derivative, 0 is the absolute orientation
	 * @param value standard deviation per second in rad/s, rad/s^2, rad/s^3, ..
	 */
	void setOriPN( const std::size_t order, const double value )
	{
		if ( order == 0 )
			for ( std::size_t i ( 3 + 3 * m_posOrder  ); i < 7 + 3 * m_posOrder; i++ )
				m_processNoise( i ) = value * value;
		else
			for ( std::size_t i( 4 + 3 * ( m_posOrder + order ) ); i < 7 + 3 * ( m_posOrder + order ); i++ )
				m_processNoise( i ) = value * value;
	}

	/**
	 * applies the process noise to a covariance matrix.
	 * @param cov the covariance matrix to change
	 * @param deltaTime time (in seconds)
	 */
	template< class M >
	void addNoise( M& cov, double deltaTime ) const
	{
		double dt = fabs( deltaTime );
		for ( int i = 0; i < stateSize(); i++ )
			cov( i, i ) += dt * m_processNoise( i );
	}
	
protected:
	int m_posOrder;
	int m_oriOrder;
	Math::Vector< double > m_processNoise;
};

} } // namespace Ubitrack::Tracking

#endif

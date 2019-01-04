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
* Class for kalman filtering of poses
*
* @author Frieder Pankratz <pankratz@in.tum.de>
*/

#ifndef __UBITRACK_TRACKING_OUTORORDERPOSEKALMANFILTER_H_INCLUDED__
#define __UBITRACK_TRACKING_OUTORORDERPOSEKALMANFILTER_H_INCLUDED__


#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/ErrorVector.h>
#include <utMeasurement/Measurement.h>
#include <utTracking/PoseKalmanFilter.h>
#include <boost/scoped_ptr.hpp>
#include <log4cpp/Category.hh>
#include <boost/numeric/ublas/io.hpp>


#include <sstream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/circular_buffer.hpp>

namespace Ubitrack {
	namespace Tracking {

		
		class UBITRACK_EXPORT OutOfOrderPoseKalmanFilter
		{
		public:
			/**
			* Constructor.
			* @param motionModel Motion model that defines the process noise and how many derivatives for position
			* and orientation to use.
			* @param bInsideOut if true, a motion model is used that assumes a correlation between orientation and
			* translation, e.g. when a non-moveable object is tracked by a mobile camera.
			*/
			OutOfOrderPoseKalmanFilter(int historyCount, std::vector< double > posPN, std::vector< double > oriPN, bool bInsideOut = false);
			
			~OutOfOrderPoseKalmanFilter(){
				m_pKF.reset();
			}
			/**
			* integrate an absolute pose measurement.
			* @param m the measured pose with timestamp and error
			*/
			void addPoseMeasurement(const Measurement::ErrorPose& m, int index=0);

			/**
			* compute a rotation for a given time, which may lie in the future
			*/
			Measurement::ErrorPose predictPose(Measurement::Timestamp t);

			void Reset();

			/** type of internal state representation */
			typedef Math::Vector< double > StateType;

			/** type of internal state representation */
			typedef Math::Matrix< double, 0, 0 > CovarianceType;

			/** returns the internal state */
			const StateType& getState() const
			{
				return m_pKF->getState();
			}

			/** returns the internal covariance */
			const CovarianceType& getCovariance() const
			{
				return m_pKF->getCovariance();
			}

			/** returns the motion model */
			const LinearPoseMotionModel& getMotionModel() const
			{
				return m_pKF->getMotionModel();
			}


		protected:

			// the kalman filter
			boost::scoped_ptr< Tracking::PoseKalmanFilter > m_pKF;


			///** output stream */
			//boost::scoped_ptr< std::ofstream > m_pStream;

			//boost::scoped_ptr< std::ofstream > m_pStreamPose;
			///** output archive */
			//boost::scoped_ptr< boost::archive::text_oarchive > m_pArchive;

			//boost::scoped_ptr< std::ofstream > m_pStreamPoseInput;
			//boost::scoped_ptr< boost::archive::text_oarchive > m_pArchiveInput;

			boost::circular_buffer<Measurement::ErrorPose> m_history;
			std::vector< double > m_posPN;
			std::vector< double > m_oriPN;
			bool m_bInsideOut;
		};

	}
} // namespace Ubitrack::Tracking

#endif // HAVE_LAPACK

#endif

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
* Component for kalman filtering
*
* @author Frieder Pankratz <pankratz@in.tum.de>
*/

#include "OutOfOrderPoseKalmanFilter.h"
#include <utUtil/Exception.h>
// get a logger
static log4cpp::Category& logger(log4cpp::Category::getInstance("Ubitrack.Tracking.OutOfOrderKalmanFilter"));

namespace Ubitrack {
	namespace Tracking {
		OutOfOrderPoseKalmanFilter::OutOfOrderPoseKalmanFilter(int historyCount, std::vector< double > posPN, std::vector< double > oriPN, bool bInsideOut){
			m_posPN = posPN;
			m_oriPN = oriPN;
			m_bInsideOut = bInsideOut;

			m_history = boost::circular_buffer<Measurement::ErrorPose>(historyCount);

			Reset();


			//// create ofstream
			//std::string sFilename = "KalmanLog.txt";

			//

			//m_pStream.reset(new std::ofstream(sFilename.c_str()));
			//if (!m_pStream->good())
			//	UBITRACK_THROW("Could not open file " + sFilename + " for writing");

			//// create oarchive
			//std::string sFilenameOut = "KalmanOutputLog";
			//for (int i = 0; i < posPN.size(); i++) {
			//	sFilenameOut += "_" + std::to_string(posPN[i]);
			//}
			//sFilenameOut += ".txt";

			//m_pStreamPose.reset(new std::ofstream(sFilenameOut.c_str()));
			//if (!m_pStreamPose->good())
			//	UBITRACK_THROW("Could not open file " + sFilenameOut + " for writing");
			//m_pArchive.reset(new boost::archive::text_oarchive(*m_pStreamPose));


			//std::string sFilenameInOut = "KalmanInputLog.txt";
			//m_pStreamPoseInput.reset(new std::ofstream(sFilenameInOut.c_str()));
			//if (!m_pStreamPoseInput->good())
			//	UBITRACK_THROW("Could not open file " + sFilenameInOut + " for writing");
			//m_pArchiveInput.reset(new boost::archive::text_oarchive(*m_pStreamPoseInput));

		}

		/**
		* integrate an absolute pose measurement.
		* @param m the measured pose with timestamp and error
		*/
		void OutOfOrderPoseKalmanFilter::addPoseMeasurement(const Measurement::ErrorPose& m, const int index){
			Measurement::Timestamp tsNow = Measurement::now();
			/*std::string linesep("\n");
			(*m_pStream) << linesep;
			(*m_pStream) << index;
			(*m_pStream) << " ";
			(*m_pStream) << tsNow;
			(*m_pStream) << " ";
			(*m_pStream) << m.time();

			(*m_pArchiveInput) << linesep;
			(*m_pArchiveInput) << m;*/

			
			bool outOfOrder = false;

			//const double stdDevPos = std::sqrt(m->covariance()(0, 0) + m->covariance()(1, 1) + m->covariance()(2, 2));

			//if (stdDevPos > 0.1) {
			//	LOG4CPP_WARN(logger, "std dev too big skip measurement, index : stdDev  " << index << " : " << stdDevPos);
			//	return;
			//}

			if (m_history.size() == 0) {
				m_pKF->addPoseMeasurement(m);
				m_history.push_back(m);
				//LOG4CPP_WARN(logger, "in m_history " << m_history.size() << " : " << m.time());
				return;
			}

			Measurement::Timestamp timeBack = m_history.rbegin()->time();
			
			int maxTimeDiffForSync = 8;

			int timeDiff = 0;

			if (m.time() > timeBack)
				timeDiff = (m.time() - timeBack) / 1000000l;
			else {
				timeDiff = (timeBack - m.time()) / 1000000l;
				timeDiff = -timeDiff;
			}
				
			
			//Measurement::Timestamp timediff = m_history.rbegin()->time() - m.time();

			//LOG4CPP_WARN(logger, "timeDiff " << timeDiff << " : " << m.time());
			// is new measurment after last input?
			if (m.time() > timeBack && timeDiff > maxTimeDiffForSync) {
				// everything ok
				LOG4CPP_WARN(logger, "new measurement index:timeDiff " << index << " : " << timeDiff << " : " << m.time() / 1000000l);
				
				m_pKF->addPoseMeasurement(m);
				m_history.push_back(m);
				return;
			}

			if (std::abs(timeDiff) <= maxTimeDiffForSync) {
				LOG4CPP_WARN(logger, "timediff small asume sync index:timeDiff " << index << " : " << timeDiff << " : " << m.time() / 1000000l);
				Measurement::ErrorPose newM(timeBack, *m);
				m_pKF->addPoseMeasurement(newM);
				m_history.push_back(newM);
				return;
			}

			// resset kalman filter

			//LOG4CPP_WARN(logger, "out of order " << m_history.back().time() << " : " << m.time());

		/*	for (int i = m_history.size() - 1; i >= 0; --i) {
				LOG4CPP_WARN(logger, "m_history " << m_history[i].time());
			}
			LOG4CPP_WARN(logger, "-------------------------------------");
			for (boost::circular_buffer<Measurement::ErrorPose>::iterator it = m_history.begin(); it != m_history.end(); ++it) {
				LOG4CPP_WARN(logger, "m_history " << it->time());
				
			}*/

			int countIndex = 0;
			bool posFound = false;
			


			if (m.time() < m_history.begin()->time()){
				LOG4CPP_WARN(logger, "event too old, reject: " << timeDiff);
				return;
			}
			
			
			if (std::abs(timeDiff) > 100){
				LOG4CPP_WARN(logger, "timdiff too big index:timediff " << index << " : " << timeDiff);
				return;
			}


			for (boost::circular_buffer<Measurement::ErrorPose>::iterator it = m_history.begin(); it != m_history.end(); ++it) {
				
				if (m.time() < it->time()) {
					
					LOG4CPP_WARN(logger, "position found for " << m.time() << " : " << countIndex << " bufferSize: " << m_history.size() << " diff to newest: " << timeDiff);

					m_history.insert(it, m);
					posFound = true;
					if (m_history.size() > 3){
						//LOG4CPP_WARN(logger, "now order like this: " << m_history[countIndex - 1].time()  << " : " << m_history[countIndex].time() << " : " << m_history[countIndex + 1].time());
					}
					
					break;
				}
				countIndex++;
			}

			if (posFound)
				Reset();
								
				

		}

		Measurement::ErrorPose OutOfOrderPoseKalmanFilter::predictPose(Measurement::Timestamp t){
			Measurement::Timestamp tsNow = Measurement::now();
			LOG4CPP_DEBUG(logger, "Computing pose for t=" << t);
			LOG4CPP_TRACE(logger, "state: " << m_pKF->getState() << std::endl << m_pKF->getCovariance());

			Measurement::ErrorPose result = m_pKF->predictPose(t);

			/*std::string linesep("\n");
			(*m_pStream) << linesep;
			(*m_pStream) << -1;
			(*m_pStream) << " ";
			(*m_pStream) << tsNow;
			(*m_pStream) << " ";
			(*m_pStream) << t;


			(*m_pArchive) << linesep;
			(*m_pArchive) << result;*/

			return result;
		}


		void OutOfOrderPoseKalmanFilter::Reset() {
			Tracking::LinearPoseMotionModel motionModel(m_posPN.size() - 1, m_oriPN.size() - 1);
			for (std::size_t i(0); i < m_posPN.size(); i++)
				motionModel.setPosPN(i, m_posPN[i]);
			for (std::size_t i(0); i < m_oriPN.size(); i++)
				motionModel.setOriPN(i, m_oriPN[i]);

			// initialize kalman filter with motion model
			m_pKF.reset(new Tracking::PoseKalmanFilter(motionModel, m_bInsideOut));

			for (boost::circular_buffer<Measurement::ErrorPose>::iterator it = m_history.begin(); it != m_history.end(); ++it) {
				m_pKF->addPoseMeasurement(*it);

			}

		}

		
	}
} // namespace Ubitrack::Tracking

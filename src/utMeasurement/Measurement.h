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
 * @ingroup datastructures
 * @file
 * Measurement template class and several default instanciations
 * @author Florian Echtler <echtler@in.tum.de>
 */


#ifndef _Ubitrack_Measurement_Measurement_INCLUDED_
#define _Ubitrack_Measurement_Measurement_INCLUDED_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/shared_ptr.hpp>

#include <iostream>
#include <vector>

#include "Timestamp.h"

#include <utMath/Vector.h>
#include <utMath/Quaternion.h>
#include <utMath/Matrix.h>
#include <utMath/Pose.h>
#include <utMath/ErrorPose.h>
#include <utMath/Scalar.h>
#include <utMath/RotationVelocity.h>
#include <utMath/CameraIntrinsics.h>

namespace Ubitrack { namespace Measurement {

// forward declaration of Measurement
template< typename Type > class Measurement;

//single measurements
typedef Measurement< Math::Scalar< double > > Distance;
typedef Measurement< Math::Vector< 2 >      > Position2D;
typedef Measurement< Math::Vector< 3 >      > Position;
typedef Measurement< Math::Vector< 3 >      > Vector3D; // e.g. magnetic field vector
typedef Measurement< Math::Vector< 4 >      > Vector4D;
typedef Measurement< Math::Vector< 8 >      > Vector8D; // New distortion model
typedef Measurement< Math::Quaternion       > Rotation;
typedef Measurement< Math::Matrix< 3, 3 >   > Matrix3x3;
typedef Measurement< Math::Matrix< 3, 4 >   > Matrix3x4;
typedef Measurement< Math::Matrix< 4, 4 >   > Matrix4x4;
typedef Measurement< Math::Scalar< int >    > Button;
typedef Measurement< Math::Pose             > Pose;
typedef Measurement< Math::ErrorPose        > ErrorPose;
typedef Measurement< Math::ErrorVector< 3 > > ErrorPosition;
typedef Measurement< Math::RotationVelocity > RotationVelocity;
typedef Measurement< Math::CameraIntrinsics< double > > CameraIntrinsics;

//multiple measurements
typedef Measurement< std::vector < Math::Pose > > PoseList;
typedef Measurement< std::vector < Math::Vector< 2 > > > PositionList2;
typedef Measurement< std::vector < Math::Vector< 3 > > > PositionList;
typedef Measurement< std::vector < Math::Scalar< double > > > DistanceList;
typedef Measurement< std::vector < Math::Scalar< unsigned long > > > IDList;
typedef Measurement< std::vector < Math::ErrorPose > > ErrorPoseList;
typedef Measurement< std::vector < Math::ErrorVector< 2 > > > ErrorPositionList2;
typedef Measurement< std::vector < Math::ErrorVector< 3 > > > ErrorPositionList;

//typedef Measurement< Math::ErrorFeaturePosition< 3 > > ErrorFeaturePosition;
//typedef Measurement< std::vector < Math::ErrorFeaturePosition< 3 , FeatureDescriptor> > > ErrorFeaturePositionList3D<class FeatureDescriptor> ;
//typedef Measurement< std::vector < Math::ErrorFeaturePosition< 2 , class FeatureDescriptor> > > ErrorFeaturePositionList2D;

/**
 * stream output operator
 */
template< typename Type > 
std::ostream& operator<< ( std::ostream& s, const Measurement< Type >& m )
{
    if ( m.invalid() ) {
        return s << "INVALID";
    } else {
	    return s << *m << " " << timestampToShortString( m.m_timestamp );
    }
}


/**
 * @ingroup datastructures
 * Measurement: template class for all measurements, usually derived from Math type.
 * M. aggregates the measurement value and a timestamp, provides (un-)serialization
 * and human-readable output.
 *
 * The \c Measurement class is derived from \c boost::shared_ptr< Type >, so use 
 * pointer syntax to access the payload. Also notice that copying the measurement
 * does not copy data, but creates a reference to the old measurement, so
 * @verbatim 
Measurement< X > a( t, dataD ); 
Measurement< X > b = a; 
*a = dataE; 
@endverbatim
 * results in \c b being changed to \c dataE, too! Use \c Measurement::clone 
 * to copy data instead.
 *
 * @param Type data type of payload.
 */
template< typename Type > 
class Measurement
	: public boost::shared_ptr< Type >
{
	friend std::ostream& operator<< <> ( std::ostream& s, const Measurement< Type >& m );
	friend class ::boost::serialization::access;

	public:
		/** remembers the content type for templated functions and classes */
		typedef Type value_type;

		/**
		 * Default Constructor
         * At least, init timestamp to zero and hence flag as invalid
		 */
		Measurement( )
            : m_timestamp ( INVALID )
		{ }

		/** Construct from timestamp. The payload is empty. */
		explicit Measurement( Timestamp t )
			: m_timestamp( t )
		{ }

		/** Construct from payload \c shared_ptr, with timestamp of 0. */
		explicit Measurement( boost::shared_ptr< Type > p )
			: boost::shared_ptr< Type>( p )
			, m_timestamp( 0 )
		{ }
		
		/** Construct from payload reference (content will be copied), with timestamp of 0. */
		explicit Measurement( const Type& m )
			: boost::shared_ptr< Type>( new Type( m ) )
			, m_timestamp( 0 )
		{ }

		/** Construct from timestamp and payload \c shared_ptr. */
		Measurement( Timestamp t, boost::shared_ptr< Type > p )
			: boost::shared_ptr< Type>( p )
			, m_timestamp( t )
		{ }

		/** Construct from timestamp and payload reference (content will be copied). */
		Measurement( Timestamp t, const Type& m )
			: boost::shared_ptr< Type>( new Type( m ) )
			, m_timestamp( t )
		{ }

		/**
		 * set the internal timestamp
		 */
		void time( Timestamp t )
		{ m_timestamp = t; }

		/**
		 * get the internal timestamp
		 */
		Timestamp time() const
		{ return m_timestamp; }
		
		/**
		 * returns a measurement that does NOT hold a reference to this object
		 */
		Measurement clone() const
		{ if ( this->get() != 0 )
            return Measurement( m_timestamp, *(this->get()) ); 
        else
            return Measurement( m_timestamp );
        }

        /**
         * Checks, if measurement is valid
         */
        bool invalid( void ) const
        {
            return m_timestamp == INVALID;
        }

        /**
         * Sets the current measurement as invalid
         */
        void invalidate( void ) {
            m_timestamp = INVALID;
        }

        template< typename Type2 >
        Measurement< Type2 > constCast() {
            return Measurement< Type2 > ( m_timestamp, boost::const_pointer_cast< Type2 >( *this ) );
        }

	protected:

		/**
		 * (un-)serialization helper function
		 */
		template< class Archive > 
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & m_timestamp;
			ar & ( *(this->get()) );
		}

		Timestamp m_timestamp;

        // Define 0 as invalid measurement
        static const Timestamp INVALID = 0;
};


} } // namespace Ubitrack::Measurement



#endif // _Ubitrack_Measurement_Measurement_INCLUDED_


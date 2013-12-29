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
 * @ingroup calibration
 * @file
 * 2D-3D pose optimization component for multiple-camera systems.
 * Requires an initial pose!
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */


#include "BundleAdjustment.h"

#include <utMath/Vector.h>
#include <utUtil/Exception.h>

//#define OPTIMIZATION_LOGGING // use define before optimization functions
#include <utMath/Optimization/LevenbergMarquardt.h>

#include <memory> // std::allocator
#include <numeric> //std::accumulate
#include <iterator> // std::iterator_traits
#include <functional> // std::mem_fun_ref

// get a logger
#include <log4cpp/Category.hh>
static log4cpp::Category& logger( log4cpp::Category::getInstance( "Ubitrack.Calibration.BundleAdjustment" ) );
static log4cpp::Category& optLogger( log4cpp::Category::getInstance( "Ubitrack.Calibration.BundleAdjustment.LM" ) );


namespace Ubitrack { namespace Calibration {

// In future, this file should not be compiled at all if lapack is not available
#ifdef HAVE_LAPACK


template< class VType >
class MinimizeReprojectionErrorAllPoints
{
protected:
	const std::size_t n_cams;
	const std::size_t n_pts3D;	

public:
	MinimizeReprojectionErrorAllPoints(
		  const std::size_t cams
		, const std::size_t points
		)
		: n_cams ( cams )
		, n_pts3D ( points )
	{}

	// /**
	 // * return the size of the result vector
	 // */
	// std::size_t size() const
	// { return n_pts3D; }
	// // return only visible points
	// // { return 2 * m_vis.size(); }


	/**
	 * @param result vector to store the result in
	 * @param input containing the parameters (target pose as 7-vector)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
	
		namespace ublas = boost::numeric::ublas;
		
		LOG4CPP_DEBUG( optLogger, "Input    Vector dimension: " << input.size() << "\n" );
		LOG4CPP_DEBUG( optLogger, "Result   Vector dimension: " << result.size() << "\n" );
		LOG4CPP_DEBUG( optLogger, "Jacobian Matrix dimension: " << J.size1() << "x" << J.size2() << " (rows x columns)\n" );
		
		LOG4CPP_TRACE( optLogger, "Input vector: \n " << static_cast< Math::Vector< VType > >( input ) << "\n" );
		
		J = Math::Matrix< VType >::zeros( J.size1(), J.size2() );
				
		std::size_t row_index = 0;
		std::size_t index_cam = 0;
		std::size_t start_index_3d_pts = n_cams*7;
		for( std::size_t iter_c( 0 ); iter_c < n_cams; ++iter_c )// für alle cameras
		{
			const std::size_t camIndex = iter_c*7;
			const VType qx = input( camIndex + 0 );
			const VType qy = input( camIndex + 1 );
			const VType qz = input( camIndex + 2 );
			const VType qw = input( camIndex + 3 );
			const VType tx = input( camIndex + 4 );
			const VType ty = input( camIndex + 5 );
			const VType tz = input( camIndex + 6 );
			
			const Math::Pose pose( Math::Quaternion( qx, qy, qz, qw ), Math::Vector< VType, 3 >( tx, ty, tz ) );

			for ( std::size_t iter_p( 0 ); iter_p < n_pts3D; ++iter_p, row_index += 2 )
			{
				// fetch x, y and z coordinate of 3D point from vector
				const std::size_t pointIndex = start_index_3d_pts + (iter_p*3);
				const VType x = input( pointIndex + 0 );
				const VType y = input( pointIndex + 1 );
				const VType z = input( pointIndex + 2 );
				
				Math::Vector< VType, 3 > pts ( x, y, z );
				pts = pose * pts;
				result( row_index + 0 ) = pts( 0 ) / pts( 2 );
				result( row_index + 1 ) = pts( 1 ) / pts( 2 );
				
				const std::size_t jrow0 = row_index;
				const std::size_t jrow1 = jrow0 + 1;
				
				const VType t2 = qw*qw;
				const VType t3 = qx*qx;
				const VType t4 = qy*qy;
				const VType t5 = qz*qz;
				const VType t6 = qw*qy*2;
				const VType t7 = t2-t3-t4+t5;
				const VType t8 = t7*z;
				const VType t9 = qx*qz*2;
				const VType t10 = qw*qx*2;
				const VType t11 = qy*qz*2;
				const VType t12 = t10+t11;
				const VType t13 = t12*y;
				const VType t14 = t6-t9;
				const VType t23 = t14*x;
				const VType t15 = t8+t13-t23+tz;
				const VType t16 = t2+t3-t4-t5;
				const VType t17 = t16*x;
				const VType t18 = qw*qz*2;
				const VType t33 = qx*qy*2;
				const VType t19 = t18-t33;
				const VType t20 = t6+t9;
				const VType t21 = t20*z;
				const VType t34 = t19*y;
				const VType t22 = t17+t21-t34+tx;
				const VType t24 = 1/(t15*t15);
				const VType t25 = qz*x*2;
				const VType t26 = qw*y*2;
				const VType t43 = qx*z*2;
				const VType t27 = t25+t26-t43;
				const VType t28 = 1/t15;
				const VType t29 = qx*x*2;
				const VType t30 = qy*y*2;
				const VType t31 = qz*z*2;
				const VType t32 = t29+t30+t31;
				const VType t35 = qw*x*2;
				const VType t36 = qy*z*2;
				const VType t44 = qz*y*2;
				const VType t37 = t35+t36-t44;
				const VType t38 = qx*y*2;
				const VType t39 = qw*z*2;
				const VType t41 = qy*x*2;
				const VType t40 = t38+t39-t41;
				const VType t42 = t28*t40;
				const VType t45 = t2-t3+t4-t5;
				const VType t46 = t45*y;
				const VType t47 = t18+t33;
				const VType t48 = t47*x;
				const VType t49 = t10-t11;
				const VType t52 = t49*z;
				const VType t50 = t46+t48-t52+ty;
				const VType t51 = t28*t37;
				J( jrow0, camIndex + 0 ) = t32/(t8+t13+tz-x*(t6-qx*qz*2))-t22*t24*t27;
				J( jrow0, camIndex + 1 ) = t42+t22*t24*t37;
				J( jrow0, camIndex + 2 ) = -t27*t28-t22*t24*t32;
				J( jrow0, camIndex + 3 ) = t51-t22*t24*t40;
				J( jrow0, camIndex + 4 ) = t28;
				J( jrow0, camIndex + 6 ) = -t22*t24;
				J( jrow0,  pointIndex + 0 ) = t16*t28+t14*t22*t24;
				J( jrow0,  pointIndex + 1 ) = -t19*t28-t12*t22*t24;
				J( jrow0,  pointIndex + 2 ) = t20*t28-t7*t22*t24;
				J( jrow1, camIndex + 0 ) = -t42-t24*t27*t50;
				J( jrow1, camIndex + 1 ) = t28*t32+t24*t37*t50;
				J( jrow1, camIndex + 2 ) = t51-t24*t32*t50;
				J( jrow1, camIndex + 3 ) = t27*t28-t24*t40*t50;
				J( jrow1, camIndex + 5 ) = t28;
				J( jrow1, camIndex + 6 ) = -t24*t50;
				J( jrow1,  pointIndex + 0 ) = t28*t47+t14*t24*t50;
				J( jrow1,  pointIndex + 1 ) = t28*t45-t12*t24*t50;
				J( jrow1,  pointIndex + 2 ) = -t28*t49-t7*t24*t50;
			}
		}
		
		LOG4CPP_TRACE( optLogger, "Result " << static_cast< Math::Vector< VType > > ( result ) << "\n" );
		LOG4CPP_TRACE( optLogger, "Jacobian " << static_cast< Math::Matrix< VType > > ( J ) << "\n" );
	}
};

/** 
 * @tparam ForwardIterator1 iterator to container including containers of 2D observations
 * @tparam ForwardIterator2 iterator to container including extrinsic camera pose
 * @tparam ForwardIterator3 iterator to container including 3D points
 */
template< typename ForwardIterator1, typename ForwardIterator2, typename ForwardIterator3 >
void simpleBundleAdjustmentImpl (
	  const ForwardIterator1 i2DPtsBegin // e.g. std::vector< std::vector< Math::Vector< T, 2 > > >::iterator -> begin()
	, const ForwardIterator1 i2DPtsEnd  // e.g. std::vector< std::vector< Math::Vector< T, 2 > > >::iterator -> end()
	, ForwardIterator2 iExtrinsicsPose // e.g. std::vector < Math::Pose >::iterator -> begin()
	, ForwardIterator3 i3DPtsBegin //  e.g. std::vector < Math::Vector< T, 3 > >::iterator -> begin()
	, ForwardIterator3 i3DPtsEnd //  e.g. std::vector < Math::Vector< T, 3 > >::iterator -> end()
	//, visibility <- next to come :)
	)
{
	typedef typename std::iterator_traits< ForwardIterator3 >::value_type vector3d_type;
	typedef typename std::iterator_traits< ForwardIterator1 >::value_type container2d_type;
	typedef typename container2d_type::iterator container_pts_iterator;
	typedef typename container2d_type::value_type vector2d_type;
	typedef typename vector2d_type::value_type value_type;
	
	const std::size_t n_cams( std::distance( i2DPtsBegin, i2DPtsEnd ) ); // <- number of cameras should equal number of containers for observations
	const std::size_t n_pts3D( std::distance( i3DPtsBegin, i3DPtsEnd ) ); // <- can be different from number of observations
	
	LOG4CPP_DEBUG( logger, "Started BundleAdjustment with " <<  n_cams << " cameras and " << n_pts3D << " points to optimize." );
	
	// count the observations for each camera individually
	std::vector< std::size_t > point_count;
	point_count.reserve( n_cams );
	std::transform( i2DPtsBegin, i2DPtsEnd, std::back_inserter( point_count ), std::mem_fun_ref( &container2d_type::size )  );
	
	// count all observations
	const std::size_t observationCountTotal = std::accumulate( point_count.begin(), point_count.end(), 0 );
	
	LOG4CPP_DEBUG( logger, "Counted " << observationCountTotal << " observations from all "  << n_cams << " cameras, creating the measurement vector" );
	
	
	// Now create the measurement vector from the 2D observations for LM optimization
	Math::Vector< value_type > observationVector( 2 * observationCountTotal );
	ForwardIterator1 iterCams = i2DPtsBegin;
	std::size_t iIndex( 0 );
	for ( std::size_t cameraIndex = 0; cameraIndex < n_cams; ++cameraIndex, ++iterCams )
	{
		const std::size_t n_points = point_count[ cameraIndex ];
		for ( std::size_t pointIndex = 0; pointIndex < n_points; ++pointIndex, ++iIndex )
		{
			const vector2d_type pt2D = *( iterCams->begin() + pointIndex );
			boost::numeric::ublas::subrange( observationVector, 2 * iIndex, 2 * (iIndex+1) ) = pt2D;
			// LOG4CPP_TRACE( logger, "2D point #" << pointIndex << " in camera #" << cameraIndex << " : " << pt2D );
		}
	}
	
	LOG4CPP_DEBUG( logger, "size of observation vector " << observationVector.size() );
	LOG4CPP_TRACE( logger, "observation vector:\n" << observationVector );
	
		
	//create parameter vector to be optimized in the Levenberg-Marquardt Optimization:
	// there are 4 values (quaternion) and 3 values (translation for each camera) and 3 values for each 3D point:
	const std::size_t vetorSize = n_cams *( 4 + 3 ) + n_pts3D * 3;
	Math::Vector< value_type > paramVector( vetorSize );
	
	ForwardIterator2 poseIter = iExtrinsicsPose;
	for ( std::size_t cameraIndex = 0; cameraIndex < n_cams; ++cameraIndex, ++poseIter )
	{
		const std::size_t index = 7*cameraIndex;
		// boost::numeric::ublas::subrange( paramVector, index, index+4 ) = poseIter->rotation();
		Math::Vector< value_type, 4 > quatVec;
		poseIter->rotation().toVector( quatVec );
		boost::numeric::ublas::subrange( paramVector, index, index+4 ) = quatVec;
		boost::numeric::ublas::subrange( paramVector, index+4, index+7 ) = poseIter->translation();
		// LOG4CPP_TRACE( logger, "Camera #" << cameraIndex << " translation: " << poseIter->translation() << ", quaternion: " << poseIter->rotation() << "." );
	}
	
	// will be done more nicely later, promise
	ForwardIterator3 pts3D = i3DPtsBegin;
	for ( std::size_t pointIndex = 0; pointIndex < n_pts3D; ++pointIndex, ++poseIter, ++pts3D )
	{
		const std::size_t index = (n_cams*7)+3*pointIndex;
		
		boost::numeric::ublas::subrange( paramVector, index, index+3 ) = *pts3D;
		LOG4CPP_TRACE( logger, "3D point #" << pointIndex << ": " << *pts3D );
	}	
	
	LOG4CPP_DEBUG( logger, "size of parameter vector " << paramVector.size() );
	LOG4CPP_TRACE( logger, "parameter vector:\n" << paramVector );
	
	
	OPT_LOG_DEBUG( "Optimizing pose over " << numberCameras << " cameras using " << observationCountTotal << " observations" );
	MinimizeReprojectionErrorAllPoints< value_type > minimizeFunc( n_cams, n_pts3D );
	value_type res = Math::Optimization::levenbergMarquardt( minimizeFunc, paramVector, observationVector, Math::Optimization::OptTerminate( 10, 1e-6 ), Math::Optimization::OptNoNormalize() );
	
	// LOG4CPP_TRACE( logger, "optimized parameter vector:\n" << paramVector );
	
	
	pts3D = i3DPtsBegin;
	for ( std::size_t pointIndex = 0; pointIndex < n_pts3D; ++pointIndex, ++pts3D )
	{
		const std::size_t index = (n_cams*7)+pointIndex*3;
		*pts3D = boost::numeric::ublas::subrange( paramVector, index, index+3 );
		LOG4CPP_TRACE( logger, "3D point #" << pointIndex << ": " << *pts3D << " (updated)" );
	}
	
	poseIter = iExtrinsicsPose;
	for ( std::size_t cameraIndex = 0; cameraIndex < n_cams; ++cameraIndex, ++poseIter )
	{
		const std::size_t index = 7*cameraIndex;
		const Math::Quaternion quat = Math::Quaternion::fromVector( boost::numeric::ublas::subrange( paramVector, index, index+4 ) ).normalize();
		const Math::Vector< value_type, 3 > trans ( boost::numeric::ublas::subrange( paramVector, index+4, index+7 ) );
		*poseIter = Math::Pose( quat, trans );
		LOG4CPP_TRACE( logger, "Camera #" << cameraIndex << " translation: " << poseIter->translation() << ", quaternion: " << poseIter->rotation() << " (updated)." );
	}
};

void simpleBundleAdjustment( const std::vector< std::vector< Math::Vector2d > >& pts2D, std::vector< Math::Pose >& poses, std::vector< Math::Vector3d > & pts3D )
{
	simpleBundleAdjustmentImpl( pts2D.begin(), pts2D.end(), poses.begin(), pts3D.begin(), pts3D.end() );
}

void simpleBundleAdjustment( const std::vector< std::vector< Math::Vector2f > >& pts2D,  std::vector< Math::Pose >& poses, std::vector< Math::Vector3f > & pts3D )
{
	simpleBundleAdjustmentImpl( pts2D.begin(), pts2D.end(), poses.begin(), pts3D.begin(), pts3D.end() );
}

#endif // HAVE_LAPACK

} } // namespace Ubitrack::Calibration

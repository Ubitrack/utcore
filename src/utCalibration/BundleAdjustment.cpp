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
 * @author Christian Waechter <daniel.pustka@in.tum.de>
 * @author Daniel Pustka <daniel.pustka@in.tum.de> (basis)
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


template< class VType, typename IntrinsicsIterator >
class MinimizeReprojectionErrorAllPoints
{
protected:
	const std::size_t n_cams;
	const std::size_t n_pts3D;
	const IntrinsicsIterator iterCamMat;
	

public:
	MinimizeReprojectionErrorAllPoints(
		  const std::size_t cams
		, const std::size_t points
		, const IntrinsicsIterator iterMat
		)
		: n_cams ( cams )
		, n_pts3D ( points )
		, iterCamMat( iterMat )
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
		
		LOG4CPP_DEBUG( logger, "Input    Vector dimension: " << input.size() << "\n" );
		LOG4CPP_DEBUG( logger, "Result   Vector dimension: " << result.size() << "\n" );
		LOG4CPP_DEBUG( logger, "Jacobian Matrix dimension: " << J.size1() << "x" << J.size2() << " (rows x columns)\n" );
		
		J = Math::Matrix< VType >::zeros( J.size1(), J.size2() );
		
		// LOG4CPP_DEBUG( logger, "Input Vector " <<  static_cast< Math::Vector< VType > > ( input ) << "\n" );
		// LOG4CPP_DEBUG( logger, "Result before " << static_cast< Math::Vector< VType > > ( result ) << "\n" );
		// LOG4CPP_DEBUG( logger, "Jacobian " << static_cast< Math::Matrix< VType > > ( J ) << "\n" );
		
		std::size_t index_cam = 0;
		std::size_t start_index_3d_pts = n_cams*7;
		std::size_t index_pts2 = 0;
		IntrinsicsIterator iterCam = iterCamMat;
		// for( std::size_t i_c( 0 ); i_c < n_cams; ++i_c, ++iterCam )// für alle cameras
		// {
			
			// for ( std::size_t i_p( 0 ); i_p < n_pts3D; ++i_p, ++index_pts2 ) // projetziere alle punkte...
			// {
				// const std::size_t ip3 = start_index_3d_pts + (i_p*3);
				// Math::Vector< VType, 3 > pts = ublas::subrange( input, ip3, ip3+3 );
				// pts = pose * pts;
				// pts = ublas::prod( *iterCam, pts );
				// ublas::subrange( result, index_pts2*2, (index_pts2*2)+2 ) = Math::Vector< VType, 2 >( pts( 0 ) / pts( 2 ), pts( 1 ) / pts( 2 ) );
			// }
		// }
		// LOG4CPP_DEBUG( logger, "Result after " << static_cast< Math::Vector< VType > > ( result ) << "\n" );
		
		
		
		iterCam = iterCamMat;
		std::size_t row_index = 0;
		for( std::size_t iter_c( 0 ); iter_c < n_cams; ++iter_c, ++iterCam )// für alle cameras
		{
			const std::size_t camIndex = iter_c*7;
			const VType qx = input( camIndex + 0 );
			const VType qy = input( camIndex + 1 );
			const VType qz = input( camIndex + 2 );
			const VType qw = input( camIndex + 3 );
			const VType tx = input( camIndex + 4 );
			const VType ty = input( camIndex + 5 );
			const VType tz = input( camIndex + 6 );
			
			const VType fx = (*iterCam)( 0, 0 );
			const VType fy = (*iterCam)( 1, 1 );
			const VType cx = (*iterCam)( 0, 2 );
			const VType cy = (*iterCam)( 1, 2 );
			const VType cz = (*iterCam)( 2, 2 );
			
			const Math::Pose pose( Math::Quaternion( qx, qy, qz, qw ), Math::Vector< VType, 3 >( tx, ty, tz ) );
			
			// std::cout << "Camera Rotation " << qx << " " << qy << " " << qz << " " << qw << std::endl;
			// std::cout << "Camera Translation " << tx << " " << ty << " " << tz << std::endl;
			// std::cout << "Camera Matrix " << fx << " " << fy << " " << cx << " " << cy << std::endl;
			
			for ( std::size_t iter_p( 0 ); iter_p < n_pts3D; ++iter_p, row_index += 2 )
			{
				// fetch x, y and z coordinate of 3D point from vector
				const std::size_t pointIndex = start_index_3d_pts + (iter_p*3);
				const VType x = input( pointIndex + 0 );
				const VType y = input( pointIndex + 1 );
				const VType z = input( pointIndex + 2 );
				const VType w = 1;//input( pointIndex + 3 );
				
				Math::Vector< VType, 3 > pts ( x, y, z );
				pts = pose * pts;
				pts = ublas::prod( *iterCam, pts );
				result( row_index + 0 ) = pts( 0 ) / pts( 2 );
				result( row_index + 1 ) = pts( 1 ) / pts( 2 );
				
				
				const std::size_t jrow0 = row_index;
				const std::size_t jrow1 = jrow0 + 1;
				
				
				J( jrow0, camIndex + 0 ) = (x*(cx*qz*2.0+fx*qx*2.0)+y*(cx*qw*2.0+fx*qy*2.0)-z*(cx*qx*2.0-fx*qz*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)-(cz*qz*x*2.0+cz*qw*y*2.0-cz*qx*z*2.0)*(x*(fx*(qw*qw+qx*qx-qy*qy-qz*qz)-cx*(qw*qy*2.0-qx*qz*2.0))+z*(cx*(qw*qw-qx*qx-qy*qy+qz*qz)+fx*(qw*qy*2.0+qx*qz*2.0))+y*(cx*(qw*qx*2.0+qy*qz*2.0)-fx*(qw*qz*2.0-qx*qy*2.0))+w*(cx*tz+fx*tx))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow0, camIndex + 1 ) = -(x*(cx*qw*2.0+fx*qy*2.0)-y*(cx*qz*2.0+fx*qx*2.0)+z*(cx*qy*2.0-fx*qw*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)+(cz*qw*x*2.0-cz*qz*y*2.0+cz*qy*z*2.0)*(x*(fx*(qw*qw+qx*qx-qy*qy-qz*qz)-cx*(qw*qy*2.0-qx*qz*2.0))+z*(cx*(qw*qw-qx*qx-qy*qy+qz*qz)+fx*(qw*qy*2.0+qx*qz*2.0))+y*(cx*(qw*qx*2.0+qy*qz*2.0)-fx*(qw*qz*2.0-qx*qy*2.0))+w*(cx*tz+fx*tx))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow0, camIndex + 2 ) = (x*(cx*qx*2.0-fx*qz*2.0)+y*(cx*qy*2.0-fx*qw*2.0)+z*(cx*qz*2.0+fx*qx*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)-(cz*qx*x*2.0+cz*qy*y*2.0+cz*qz*z*2.0)*(x*(fx*(qw*qw+qx*qx-qy*qy-qz*qz)-cx*(qw*qy*2.0-qx*qz*2.0))+z*(cx*(qw*qw-qx*qx-qy*qy+qz*qz)+fx*(qw*qy*2.0+qx*qz*2.0))+y*(cx*(qw*qx*2.0+qy*qz*2.0)-fx*(qw*qz*2.0-qx*qy*2.0))+w*(cx*tz+fx*tx))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow0, camIndex + 3 ) = (-x*(cx*qy*2.0-fx*qw*2.0)+y*(cx*qx*2.0-fx*qz*2.0)+z*(cx*qw*2.0+fx*qy*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)-(cz*qy*x*-2.0+cz*qx*y*2.0+cz*qw*z*2.0)*(x*(fx*(qw*qw+qx*qx-qy*qy-qz*qz)-cx*(qw*qy*2.0-qx*qz*2.0))+z*(cx*(qw*qw-qx*qx-qy*qy+qz*qz)+fx*(qw*qy*2.0+qx*qz*2.0))+y*(cx*(qw*qx*2.0+qy*qz*2.0)-fx*(qw*qz*2.0-qx*qy*2.0))+w*(cx*tz+fx*tx))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow0, camIndex + 4 ) = (fx*w)/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w);
				J( jrow0, camIndex + 6 ) = -(fx*w*(tx*w+(qw*qw)*x+(qx*qx)*x-(qy*qy)*x-(qz*qz)*x-qw*qz*y*2.0+qx*qy*y*2.0+qw*qy*z*2.0+qx*qz*z*2.0)*1.0/pow(tz*w+(qw*qw)*z-(qx*qx)*z-(qy*qy)*z+(qz*qz)*z-qw*qy*x*2.0+qx*qz*x*2.0+qw*qx*y*2.0+qy*qz*y*2.0,2.0))/cz;
				J( jrow0, pointIndex + 0 ) = (fx*1.0/pow(tz*w+(qw*qw)*z-(qx*qx)*z-(qy*qy)*z+(qz*qz)*z-qw*qy*x*2.0+qx*qz*x*2.0+qw*qx*y*2.0+qy*qz*y*2.0,2.0)*((qw*qw*qw*qw)*z-(qx*qx*qx*qx)*z+(qy*qy*qy*qy)*z-(qz*qz*qz*qz)*z+qw*(qx*qx*qx)*y*2.0+(qw*qw*qw)*qx*y*2.0-qy*(qz*qz*qz)*y*2.0-(qy*qy*qy)*qz*y*2.0+(qw*qw)*tz*w+(qx*qx)*tz*w-(qy*qy)*tz*w-(qz*qz)*tz*w+(qw*qw)*(qy*qy)*z*2.0-(qx*qx)*(qz*qz)*z*2.0+qw*qy*tx*w*2.0-qx*qz*tx*w*2.0+qw*qx*(qy*qy)*y*2.0+qw*qx*(qz*qz)*y*2.0-(qw*qw)*qy*qz*y*2.0-(qx*qx)*qy*qz*y*2.0))/cz;
				J( jrow0, pointIndex + 1 ) = (cx*(qw*qx*2.0+qy*qz*2.0)-fx*(qw*qz*2.0-qx*qy*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)-cz*(qw*qx*2.0+qy*qz*2.0)*(x*(fx*(qw*qw+qx*qx-qy*qy-qz*qz)-cx*(qw*qy*2.0-qx*qz*2.0))+z*(cx*(qw*qw-qx*qx-qy*qy+qz*qz)+fx*(qw*qy*2.0+qx*qz*2.0))+y*(cx*(qw*qx*2.0+qy*qz*2.0)-fx*(qw*qz*2.0-qx*qy*2.0))+w*(cx*tz+fx*tx))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow0, pointIndex + 2 ) = (fx*1.0/pow(tz*w+(qw*qw)*z-(qx*qx)*z-(qy*qy)*z+(qz*qz)*z-qw*qy*x*2.0+qx*qz*x*2.0+qw*qx*y*2.0+qy*qz*y*2.0,2.0)*(-(qw*qw*qw*qw)*x+(qx*qx*qx*qx)*x-(qy*qy*qy*qy)*x+(qz*qz*qz*qz)*x+qw*(qz*qz*qz)*y*2.0+qx*(qy*qy*qy)*y*2.0+(qw*qw*qw)*qz*y*2.0+(qx*qx*qx)*qy*y*2.0-(qw*qw)*tx*w+(qx*qx)*tx*w+(qy*qy)*tx*w-(qz*qz)*tx*w-(qw*qw)*(qy*qy)*x*2.0+(qx*qx)*(qz*qz)*x*2.0+qw*qy*tz*w*2.0+qx*qz*tz*w*2.0+(qw*qw)*qx*qy*y*2.0+qw*(qx*qx)*qz*y*2.0+qw*(qy*qy)*qz*y*2.0+qx*qy*(qz*qz)*y*2.0))/cz;
				J( jrow1, camIndex + 0 ) = (x*(cy*qz*2.0+fy*qy*2.0)+y*(cy*qw*2.0-fy*qx*2.0)-z*(cy*qx*2.0+fy*qw*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)-(cz*qz*x*2.0+cz*qw*y*2.0-cz*qx*z*2.0)*(y*(fy*(qw*qw-qx*qx+qy*qy-qz*qz)+cy*(qw*qx*2.0+qy*qz*2.0))+z*(cy*(qw*qw-qx*qx-qy*qy+qz*qz)-fy*(qw*qx*2.0-qy*qz*2.0))-x*(cy*(qw*qy*2.0-qx*qz*2.0)-fy*(qw*qz*2.0+qx*qy*2.0))+w*(cy*tz+fy*ty))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow1, camIndex + 1 ) = -(x*(cy*qw*2.0-fy*qx*2.0)-y*(cy*qz*2.0+fy*qy*2.0)+z*(cy*qy*2.0-fy*qz*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)+(cz*qw*x*2.0-cz*qz*y*2.0+cz*qy*z*2.0)*(y*(fy*(qw*qw-qx*qx+qy*qy-qz*qz)+cy*(qw*qx*2.0+qy*qz*2.0))+z*(cy*(qw*qw-qx*qx-qy*qy+qz*qz)-fy*(qw*qx*2.0-qy*qz*2.0))-x*(cy*(qw*qy*2.0-qx*qz*2.0)-fy*(qw*qz*2.0+qx*qy*2.0))+w*(cy*tz+fy*ty))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow1, camIndex + 2 ) = (x*(cy*qx*2.0+fy*qw*2.0)+y*(cy*qy*2.0-fy*qz*2.0)+z*(cy*qz*2.0+fy*qy*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)-(cz*qx*x*2.0+cz*qy*y*2.0+cz*qz*z*2.0)*(y*(fy*(qw*qw-qx*qx+qy*qy-qz*qz)+cy*(qw*qx*2.0+qy*qz*2.0))+z*(cy*(qw*qw-qx*qx-qy*qy+qz*qz)-fy*(qw*qx*2.0-qy*qz*2.0))-x*(cy*(qw*qy*2.0-qx*qz*2.0)-fy*(qw*qz*2.0+qx*qy*2.0))+w*(cy*tz+fy*ty))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow1, camIndex + 3 ) = (-x*(cy*qy*2.0-fy*qz*2.0)+y*(cy*qx*2.0+fy*qw*2.0)+z*(cy*qw*2.0-fy*qx*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)-(cz*qy*x*-2.0+cz*qx*y*2.0+cz*qw*z*2.0)*(y*(fy*(qw*qw-qx*qx+qy*qy-qz*qz)+cy*(qw*qx*2.0+qy*qz*2.0))+z*(cy*(qw*qw-qx*qx-qy*qy+qz*qz)-fy*(qw*qx*2.0-qy*qz*2.0))-x*(cy*(qw*qy*2.0-qx*qz*2.0)-fy*(qw*qz*2.0+qx*qy*2.0))+w*(cy*tz+fy*ty))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow1, camIndex + 5 ) = (fy*w)/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w);
				J( jrow1, camIndex + 6 ) = -(fy*w*(ty*w+(qw*qw)*y-(qx*qx)*y+(qy*qy)*y-(qz*qz)*y+qw*qz*x*2.0+qx*qy*x*2.0-qw*qx*z*2.0+qy*qz*z*2.0)*1.0/pow(tz*w+(qw*qw)*z-(qx*qx)*z-(qy*qy)*z+(qz*qz)*z-qw*qy*x*2.0+qx*qz*x*2.0+qw*qx*y*2.0+qy*qz*y*2.0,2.0))/cz;
				J( jrow1, pointIndex + 0 ) = -(cy*(qw*qy*2.0-qx*qz*2.0)-fy*(qw*qz*2.0+qx*qy*2.0))/(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w)+cz*(qw*qy*2.0-qx*qz*2.0)*(y*(fy*(qw*qw-qx*qx+qy*qy-qz*qz)+cy*(qw*qx*2.0+qy*qz*2.0))+z*(cy*(qw*qw-qx*qx-qy*qy+qz*qz)-fy*(qw*qx*2.0-qy*qz*2.0))-x*(cy*(qw*qy*2.0-qx*qz*2.0)-fy*(qw*qz*2.0+qx*qy*2.0))+w*(cy*tz+fy*ty))*1.0/pow(-cz*x*(qw*qy*2.0-qx*qz*2.0)+cz*y*(qw*qx*2.0+qy*qz*2.0)+cz*z*(qw*qw-qx*qx-qy*qy+qz*qz)+cz*tz*w,2.0);
				J( jrow1, pointIndex + 1 ) = -(fy*1.0/pow(tz*w+(qw*qw)*z-(qx*qx)*z-(qy*qy)*z+(qz*qz)*z-qw*qy*x*2.0+qx*qz*x*2.0+qw*qx*y*2.0+qy*qz*y*2.0,2.0)*(-(qw*qw*qw*qw)*z-(qx*qx*qx*qx)*z+(qy*qy*qy*qy)*z+(qz*qz*qz*qz)*z+qw*(qy*qy*qy)*x*2.0+(qw*qw*qw)*qy*x*2.0+qx*(qz*qz*qz)*x*2.0+(qx*qx*qx)*qz*x*2.0-(qw*qw)*tz*w+(qx*qx)*tz*w-(qy*qy)*tz*w+(qz*qz)*tz*w-(qw*qw)*(qx*qx)*z*2.0+(qy*qy)*(qz*qz)*z*2.0+qw*qx*ty*w*2.0+qy*qz*ty*w*2.0+qw*(qx*qx)*qy*x*2.0+(qw*qw)*qx*qz*x*2.0+qw*qy*(qz*qz)*x*2.0+qx*(qy*qy)*qz*x*2.0))/cz;
				J( jrow1, pointIndex + 2 ) = -(fy*1.0/pow(tz*w+(qw*qw)*z-(qx*qx)*z-(qy*qy)*z+(qz*qz)*z-qw*qy*x*2.0+qx*qz*x*2.0+qw*qx*y*2.0+qy*qz*y*2.0,2.0)*((qw*qw*qw*qw)*y+(qx*qx*qx*qx)*y-(qy*qy*qy*qy)*y-(qz*qz*qz*qz)*y+qw*(qz*qz*qz)*x*2.0-qx*(qy*qy*qy)*x*2.0+(qw*qw*qw)*qz*x*2.0-(qx*qx*qx)*qy*x*2.0+(qw*qw)*ty*w-(qx*qx)*ty*w-(qy*qy)*ty*w+(qz*qz)*ty*w+(qw*qw)*(qx*qx)*y*2.0-(qy*qy)*(qz*qz)*y*2.0+qw*qx*tz*w*2.0-qy*qz*tz*w*2.0-(qw*qw)*qx*qy*x*2.0+qw*(qx*qx)*qz*x*2.0+qw*(qy*qy)*qz*x*2.0-qx*qy*(qz*qz)*x*2.0))/cz;
			}
		}
		
		LOG4CPP_DEBUG( logger, "Result " << static_cast< Math::Vector< VType > > ( result ) << "\n" );
		LOG4CPP_DEBUG( logger, "Jacobian " << static_cast< Math::Matrix< VType > > ( J ) << "\n" );
	}
};

/** 
 * @tparam ForwardIterator1 iterator to container including 2D observations
 * @tparam ForwardIterator2 iterator to container including intrinsic camera matrices
 * @tparam ForwardIterator3 iterator to container including extrinsic camera pose
 * @tparam ForwardIterator4 iterator to container including 3D points
 */

template< template< typename > class ForwardIterator1, template < typename, typename > class ForwardIterator2, typename container2d_type, typename ForwardIterator3, typename ForwardIterator4, typename ForwardIterator5 >
void simpleBundleAdjustmentImpl (
	  const ForwardIterator1< typename ForwardIterator2< container2d_type, typename std::allocator< container2d_type > > > i2DPtsBegin // e.g. std::vector< std::vector< Math::Vector< T, 2 > > >::iterator -> begin()
	, const ForwardIterator1< typename ForwardIterator2< container2d_type, typename std::allocator< container2d_type > > > i2DPtsEnd  // e.g. std::vector< std::vector< Math::Vector< T, 2 > > >::iterator -> end()
	, const ForwardIterator3 iIntrinsicsMat // e.g. std::vector< Math::Matrix< T, 3, 3 > >::iterator -> begin()
	, ForwardIterator4 iExtrinsicsPose // e.g. std::vector < Math::Pose >::iterator -> begin()
	, ForwardIterator5 i3DPtsBegin //  e.g. std::vector < Math::Vector< T, 3 > >::iterator -> begin()
	, ForwardIterator5 i3DPtsEnd //  e.g. std::vector < Math::Vector< T, 3 > >::iterator -> end()
	//, visibility <- next to come :)
	)
{
	typedef typename ForwardIterator2< container2d_type, typename std::allocator< container2d_type > > container_pts2d;
	typedef typename ForwardIterator1< container_pts2d > container_pts_iterator;
	
	// typedef typename ForwardIterator1::container_type container2d_type;
	// typedef typename std::iterator_traits< ForwardIterator2 >::value_type vector2d_type;
	typedef typename container2d_type::value_type vector2d_type;
	typedef typename std::iterator_traits< ForwardIterator5 >::value_type vector3d_type;
	typedef typename vector2d_type::value_type value_type;
	
	const std::size_t n_cams( std::distance( i2DPtsBegin, i2DPtsEnd ) ); // <- number of cameras should equal number of containers for observations
	const std::size_t n_pts3D( std::distance( i3DPtsBegin, i3DPtsEnd ) ); // <- can be different from number of observations
	
	LOG4CPP_DEBUG( logger, "Started BundleAdjustment with " <<  n_cams << " cameras and " << n_pts3D << " points to optimize." );
	
	// count the observations for each camera individually
	std::vector< std::size_t > point_count;
	point_count.reserve( n_cams );
	std::transform( i2DPtsBegin, i2DPtsEnd, std::back_inserter( point_count ), std::mem_fun_ref( &container_pts_iterator::value_type::size )  );
	
	// count all observations
	const std::size_t observationCountTotal = std::accumulate( point_count.begin(), point_count.end(), 0 );
	
	LOG4CPP_DEBUG( logger, "Counted " << observationCountTotal << " observations from all "  << n_cams << " cameras, creating the measurement vector" );
	
	
	// Now create the measurement vector from the 2D observations for LM optimization
	Math::Vector< value_type > observationVector( 2 * observationCountTotal );
	container_pts_iterator iterCams = i2DPtsBegin;
	std::size_t iIndex( 0 );
	for ( std::size_t cameraIndex = 0; cameraIndex < n_cams; ++cameraIndex, ++iterCams )
	{
		const std::size_t n_points = point_count[ cameraIndex ];
		
		// typename container_pts2d::iterator iterPts = iterCams->begin() + pointIndex ;
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
	
	ForwardIterator4 poseIter = iExtrinsicsPose;
	for ( std::size_t cameraIndex = 0; cameraIndex < n_cams; ++cameraIndex, ++poseIter )
	{
		const std::size_t index = 7*cameraIndex;
		// boost::numeric::ublas::subrange( paramVector, index, index+4 ) = poseIter->rotation();
		poseIter->rotation().toVector( boost::numeric::ublas::subrange( paramVector, index, index+4 ) );
		boost::numeric::ublas::subrange( paramVector, index+4, index+7 ) = poseIter->translation();
		// LOG4CPP_TRACE( logger, "Camera #" << cameraIndex << " translation: " << poseIter->translation() << ", quaternion: " << poseIter->rotation() << "." );
	}
	
	// will be done more nicely later, promise
	ForwardIterator5 pts3D = i3DPtsBegin;
	for ( std::size_t pointIndex = 0; pointIndex < n_pts3D; ++pointIndex, ++poseIter, ++pts3D )
	{
		const std::size_t index = (n_cams*7)+3*pointIndex;
		
		boost::numeric::ublas::subrange( paramVector, index, index+3 ) = *pts3D;
		LOG4CPP_TRACE( logger, "3D point #" << pointIndex << ": " << *pts3D );
	}	
	
	LOG4CPP_DEBUG( logger, "size of parameter vector " << paramVector.size() );
	LOG4CPP_TRACE( logger, "parameter vector:\n" << paramVector );
	
	
	OPT_LOG_DEBUG( "Optimizing pose over " << numberCameras << " cameras using " << observationCountTotal << " observations" );
	MinimizeReprojectionErrorAllPoints< value_type, ForwardIterator3 > minimizeFunc( n_cams, n_pts3D, iIntrinsicsMat );
	value_type res = Math::Optimization::levenbergMarquardt( minimizeFunc, paramVector, observationVector, Math::Optimization::OptTerminate( 10, 1e-6 ), Math::Optimization::OptNoNormalize() );
	
	// LOG4CPP_TRACE( logger, "optimized parameter vector:\n" << paramVector );
	
	
	pts3D = i3DPtsBegin;
	for ( std::size_t pointIndex = 0; pointIndex < n_pts3D; ++pointIndex, ++poseIter, ++pts3D )
	{
		const std::size_t index = (n_cams*7)+3*pointIndex;
		*pts3D = boost::numeric::ublas::subrange( paramVector, index, index+3 );
		LOG4CPP_TRACE( logger, "3D point #" << pointIndex << ": " << *pts3D << " (updated)" );
	}
	
	poseIter = iExtrinsicsPose;
	for ( std::size_t cameraIndex = 0; cameraIndex < n_cams; ++cameraIndex, ++poseIter )
	{
		
		const std::size_t index = 7*cameraIndex;
		// boost::numeric::ublas::subrange( paramVector, index, index+4 ) = poseIter->rotation();
		
		const Math::Quaternion quat = Math::Quaternion::fromVector( boost::numeric::ublas::subrange( paramVector, index, index+4 ) );
		const Math::Vector< value_type, 3 > trans ( boost::numeric::ublas::subrange( paramVector, index+4, index+7 ) );
		
		*poseIter = Math::Pose( quat, trans );
		LOG4CPP_TRACE( logger, "Camera #" << cameraIndex << " translation: " << poseIter->translation() << ", quaternion: " << poseIter->rotation() << " (updated)." );
	}
	
	
};

// 
void simpleBundleAdjustment( const std::vector< std::vector< Math::Vector2d > >& pts2D, const std::vector< Math::Matrix3x3d >& mat3x3, std::vector< Math::Pose >& poses, std::vector< Math::Vector3d > & pts3D )
{
	simpleBundleAdjustmentImpl( pts2D.begin(), pts2D.end(), mat3x3.begin(), poses.begin(), pts3D.begin(), pts3D.end() );
}

void simpleBundleAdjustment( const std::vector< std::vector< Math::Vector2f > >& pts2D, const std::vector< Math::Matrix3x3f >& mat3x3, std::vector< Math::Pose >& poses, std::vector< Math::Vector3f > & pts3D )
{
	simpleBundleAdjustmentImpl( pts2D.begin(), pts2D.end(), mat3x3.begin(), poses.begin(), pts3D.begin(), pts3D.end() );
}

#endif // HAVE_LAPACK

} } // namespace Ubitrack::Calibration

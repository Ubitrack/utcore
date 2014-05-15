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
 * functions for 3D->2D projections
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include "Dehomogenization.h"
#include <utMath/Matrix.h>
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include "QuaternionRotation.h"
 
#include "RadialDistortion.h"

#ifndef __UBITRACK_CALIBRATION_FUNCTION_MULTIPLECAMERAPROJECTION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_MULTIPLECAMERAPROJECTION_H_INCLUDED__
 
namespace Ubitrack { namespace Algorithm { namespace Function {

/**
 * For a given multi-camera setup, project 3D points into each camera image plane and return 2D coordinates.
 */
template< class VType = double >
class MultipleCameraProjection
{
public:
	/** 
	 * constructor.
	 * note: all parameters must stay constant during lifetime of the object
	 * @param p3D reference to vector of 3D-points to be projected (i.e. marker positions in target coordinates
	 * @param cameraPoses reference to a vector of camera poses
	 * @param cameraIntrinsics reference to a vector of camera intrinsic parameters
	 * @param cameraDistortions reference to a vector of camera distortion parameters
	 * @param visibilities reference to a vector of observations. Each vector element contains a 
	 *    pair (i_p, i_c) which specifies that camera i_c has measured point i_p.
	 */
	MultipleCameraProjection( const std::vector< Math::Vector< VType, 3 > >& p3D, 
		const std::vector< Math::Pose >& cameraPoses, 
		const std::vector< Math::Matrix< VType, 3, 3 > >& cameraIntrinsics, 
		const std::vector< Math::Vector< VType, 4 > >& cameraDistortions, 
		const std::vector< std::pair< std::size_t, std::size_t > > visibilities )
		: m_p3D( p3D )
		, m_camP( cameraPoses )
		, m_camI( cameraIntrinsics )
		, m_camD( cameraDistortions )
		, m_vis( visibilities )
	{}

	/**
	 * return the size of the result vector
	 */
	std::size_t size() const
	{ return 2 * m_vis.size(); }

	/**
	 * @param result vector to store the result in
	 * @param input containing the parameters (target pose as 7-vector)
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		using namespace Ubitrack::Math;
		namespace ublas = boost::numeric::ublas;

		// convert quaternion to matrix (for speedup)
		Math::Quaternion rotQ( Quaternion::fromVector( ublas::subrange( input, 3, 7 ) ) );
		Math::Matrix< VType, 3, 3 > rot( rotQ );
		
		// create vectors
		Math::Vector< VType, 3 > rotated;
		Math::Vector< VType, 3 > camCoord;
		Math::Vector< VType, 2 > camCoordDehom;
		Math::Vector< VType, 2 > distorted;
		Math::Vector< VType, 2 > projected;
		
		for ( std::size_t i( 0 ); i < m_vis.size(); i++ )
		{
			// shortcuts
			const Math::Vector< VType, 3 >& p3D( m_p3D[ m_vis[ i ].first ] );
			const Math::Pose& camP( m_camP[ m_vis[ i ].second ] );
			const Math::Matrix< VType, 3, 3 >& camI( m_camI[ m_vis[ i ].second ] );
			const Math::Vector< VType, 4 >& camD( m_camD[ m_vis[ i ].second ] );
			
			// rotate & project points
			noalias( rotated ) = ublas::prod( rot, p3D ) + ublas::subrange( input, 0, 3 );
			noalias( camCoord ) = camP * rotated;
			noalias( camCoordDehom ) = ublas::subrange( camCoord, 0, 2 ) / camCoord( 2 );
			Function::radialDistortion( distorted, camCoordDehom, camD );
			noalias( projected ) = ublas::prod( ublas::subrange( camI, 0, 2, 0, 2 ), distorted );
			projected( 0 ) += camI( 0, 2 );
			projected( 1 ) += camI( 1, 2 );
			projected *= camI( 2, 2 ); // should be -1 or 1
			
			ublas::subrange( result, 2 * i, 2 * i + 2 ) = projected;
		}
	}

	/**
	 * @param result vector to store the result in
	 * @param input containing the parameters (target pose as 7-vector)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		// TODO: implement as one function (more efficient)
		evaluate( result, input );
		jacobian( input, J );
	}
	
	/**
	 * @param input containing the parameters (target pose as 7-vector)
	 * @param J matrix to store the jacobian (evaluated for param) in
	 */
	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		using namespace Ubitrack::Math;
		namespace ublas = boost::numeric::ublas;

		// convert quaternion to matrix (for speedup)
		Math::Quaternion rotQ( Quaternion::fromVector( ublas::subrange( input, 3, 7 ) ) );
		Math::Matrix< VType, 3, 3 > rot( rotQ );
		
		// create vectors
		Math::Vector< VType, 3 > rotated;
		Math::Vector< VType, 3 > camCoord;
		Math::Vector< VType, 2 > camCoordDehom;
		Math::Vector< VType, 2 > distorted;
		Math::Vector< VType, 2 > projected;
		
		// create matrices
		Math::Matrix< VType, 3, 4 > rotJ;
		Math::Matrix< VType, 2, 3 > dehomJ;
		Math::Matrix< VType, 2, 2 > distJ;
		Math::Matrix< VType, 2, 3 > projJ;

		Math::Matrix< VType, 2, 2 > jA;
		Math::Matrix< VType, 2, 3 > jB;
		Math::Matrix< VType, 2, 3 > jC;

		for ( unsigned i = 0; i < m_vis.size(); i++ )
		{
			// shortcuts
			const Math::Vector< VType, 3 >& p3D( m_p3D[ m_vis[ i ].first ] );
			const Math::Pose& camP( m_camP[ m_vis[ i ].second ] );
			const Math::Matrix< VType, 3, 3 >& camI( m_camI[ m_vis[ i ].second ] );
			const Math::Vector< VType, 4 >& camD( m_camD[ m_vis[ i ].second ] );
			
			// rotate & project points
			noalias( rotated ) = ublas::prod( rot, p3D ) + ublas::subrange( input, 0, 3 );
			noalias( camCoord ) = camP * rotated;
			noalias( camCoordDehom ) = ublas::subrange( camCoord, 0, 2 ) / camCoord( 2 );
			Function::radialDistortion( distorted, camCoordDehom, camD );
			noalias( projected ) = ublas::prod(
				ublas::subrange( camI, 0, 2, 0, 2 ), distorted );
			projected( 0 ) += camI( 0, 2 );
			projected( 1 ) += camI( 1, 2 );
			projected *= camI( 2, 2 );
			
			// create jacobian for this measurement
			QuaternionRotation< VType > qrf( p3D );
			qrf.jacobian( ublas::subrange( input, 3, 7 ), rotJ );
			Math::Matrix< VType, 3, 3 > rotCamJ( camP.rotation() );
			Dehomogenization< 3 >().jacobian( camCoord, dehomJ );
			Function::RadialDistortionWrtP< VType >( camD ).jacobian(
				camCoordDehom, distJ );

			noalias( jA ) = ublas::prod( ublas::subrange( camI, 0, 2, 0, 2 ), distJ )
				* camI( 2, 2 ); // camI( 2, 2 ) should be -1 or 1
			noalias( jB ) = ublas::prod( jA, dehomJ );
			noalias( jC ) = ublas::prod( jB, rotCamJ );
			
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 0, 3 ) ) 
				= jC;
			noalias( ublas::subrange( J, i * 2, ( i + 1 ) * 2, 3, 7 ) ) 
				= ublas::prod( jC, rotJ );
		}
	}
	
protected:
	const std::vector< Math::Vector< VType, 3 > >& m_p3D;
	const std::vector< Math::Pose >& m_camP;
	const std::vector< Math::Matrix< VType, 3, 3 > >& m_camI;
	const std::vector< Math::Vector< VType, 4 > >& m_camD;
	const std::vector< std::pair< std::size_t, std::size_t > > m_vis;
};

} } } // namespace Ubitrack::Algorithm::Function

#endif

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
 * @ingroup tracking_algorithms
 * @file
 * Implements functions for bundle adjustment.
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#include "../BundleAdjustment.h"

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <utMath/cast_assign.h>
#include "Dehomogenization.h"
#include "RadialDistortion.h"
#include "CameraIntrinsicsMultiplication.h"
#include "LieRotation.h"


namespace Ubitrack { namespace Calibration {

template< class T >
class BundleAdjustmentFunction
{
public:
	BundleAdjustmentFunction( BundleAdjustmentNetwork< T >& net )
		: m_net( net )
	{
		m_pointOffset = 0;
		m_bodyPoseOffset = m_pointOffset + 3 * m_net.points.size();
		m_imageOffset = m_bodyPoseOffset + 6 * ( m_net.bodyPoses.size() - 1 );
		m_intrinsicsOffset = m_imageOffset + 6 * m_net.images.size();
	}


	// methods part of the UnaryFunctionPrototype interface

	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& input ) const
	{
		namespace ublas = boost::numeric::ublas;
		ublas::matrix< T, ublas::column_major > J( result.size(), input.size() );
		evaluateWithJacobian( result, input, J );
	}

	template< class VT2, class MT > 
	void jacobian( const VT2& input, MT& J ) const
	{
		namespace ublas = boost::numeric::ublas;
		Math::Vector< T > result( J.size1() );
		evaluateWithJacobian( result, input, J );
	}



	/**
	 * @param result a predicted measurement vector
	 * @param input a parameter vector
	 * @param J jacobian of measurements wrt. parameters
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		namespace ublas = boost::numeric::ublas;
		// precompute image rotation matrices
		std::vector< Math::Matrix< T, 3, 3 > > camRotations( m_net.images.size() );
		std::size_t iV = m_imageOffset;
		for ( std::size_t i = 0; i != m_net.images.size(); i++, iV += 6 )
			Math::Quaternion::fromLogarithm( ublas::subrange( input, iV + 3, iV + 6 ) ).toMatrix( camRotations[ i ] );

		// precompute body rotation matrices
		std::vector< Math::Matrix< T, 3, 3 > > bodyRotations( m_net.bodyPoses.size() );
		iV = m_bodyPoseOffset;
		for ( std::size_t i = 1; i < m_net.bodyPoses.size(); i++, iV += 6 )
			Math::Quaternion::fromLogarithm( ublas::subrange( input, iV + 3, iV + 6 ) ).toMatrix( bodyRotations[ i ] );

		// clear jacobian
		noalias( J ) = ublas::zero_matrix< T >( J.size1(), J.size2() );

		std::size_t iM = 0; // offset into measurement vector

		// free point measurements
		for ( typename std::vector< typename BundleAdjustmentNetwork< T >::FreePointMeasurement >::const_iterator it = m_net.freePointMeasurements.begin();
			it != m_net.freePointMeasurements.end(); it++, iM += 2 )
		{
			std::size_t iP = m_pointOffset + 3 * it->iPoint; // offset into parameter vector

			ublas::vector_range< VT1 > resultRange( result, ublas::range( iM, iM + 2 ) );
			ublas::matrix_range< MT > jRange1( J, ublas::range( iM, iM + 2 ), ublas::range( 0, J.size2() ) );
			ublas::matrix_range< MT > jRange2( J, ublas::range( iM, iM + 2 ), ublas::range( iP, iP + 3 ) );
			Math::Vector< T, 3 > p3d( ublas::subrange( input, iP, iP + 3 ) );
			evaluateSingleWorldPointWithJacobian( resultRange, input, jRange1, jRange2, 
				camRotations, bodyRotations, it->iCamera, it->iImage, p3d );
		}

		// body point measurements
		for ( typename std::vector< typename BundleAdjustmentNetwork< T >::BodyPointMeasurement >::const_iterator it = m_net.bodyPointMeasurements.begin();
			it != m_net.bodyPointMeasurements.end(); it++, iM += 2 )
		{
			std::size_t iP = m_bodyPoseOffset + 6 * ( it->iBodyPose - 1 ); // offset into parameter vector

			// transform point from body to world
			Math::Vector< T, 3 > worldPoint;
			if ( it->iBodyPose == 0 )
			{
				// world-defining has identity-pose
				noalias( worldPoint ) = m_net.bodies[ it->iBody ][ it->iPoint ];
			}
			else
			{
				// take pose from parameter
				noalias( worldPoint ) = ublas::prod( bodyRotations[ it->iBodyPose ], m_net.bodies[ it->iBody ][ it->iPoint ] );
				noalias( worldPoint ) += ublas::subrange( input, iP, iP + 3 );
			}

			// transform from world to image and compute jacobian
			Math::Matrix< T, 2, 3 > jW2I;
			ublas::vector_range< VT1 > resultRange( result, ublas::range( iM, iM + 2 ) );
			ublas::matrix_range< MT > jRange1( J, ublas::range( iM, iM + 2 ), ublas::range( 0, J.size2() ) );
			evaluateSingleWorldPointWithJacobian( resultRange, input, jRange1, jW2I, 
				camRotations, bodyRotations, it->iCamera, it->iImage, worldPoint );

			if ( it->iBodyPose != 0 )
			{
				// jacobian of body translation
				ublas::subrange( J, iM, iM + 2, iP, iP + 3 ) = jW2I;

				// jacobian of body rotation
				Math::Matrix< T, 3, 3 > jBodyRot;
				Function::LieRotation< T >( m_net.bodies[ it->iBody ][ it->iPoint ] ).jacobian( 
					ublas::subrange( input, iP + 3, iP + 6 ), jBodyRot );
				noalias( ublas::subrange( J, iM, iM + 2, iP + 3, iP + 6 ) ) = ublas::prod( jW2I, jBodyRot );
			}
		}
	}


	// methods not part of the UnaryFunctionPrototype interface

	/** size of the measurement vector */
	std::size_t measurementSize() const
	{ return 2 * ( m_net.freePointMeasurements.size() + m_net.bodyPointMeasurements.size() ); }

	/** size of the parameter vector */
	std::size_t parameterSize() const
	{ 
		return 
			3 * m_net.points.size()
			+ 6 * ( m_net.bodyPoses.size() - 1 )
			+ 6 * m_net.images.size()
			+ ( m_net.bEstimateIntrinsics ? 9 * m_net.intrinsics.size() : 0 );
	}

	/** creates a measurement vector from the network description */
	template< class VT >
	void buildMeasurementVector( VT& v )
	{
		namespace ublas = boost::numeric::ublas;
		std::size_t iV = 0;

		// free point measurements
		for ( typename std::vector< typename BundleAdjustmentNetwork< T >::FreePointMeasurement >::iterator it = m_net.freePointMeasurements.begin();
			it != m_net.freePointMeasurements.end(); it++, iV += 2 )
			ublas::subrange( v, iV, iV + 2 ) = it->measurement;

		// body point measurements
		for ( typename std::vector< typename BundleAdjustmentNetwork< T >::BodyPointMeasurement >::iterator it = m_net.bodyPointMeasurements.begin();
			it != m_net.bodyPointMeasurements.end(); it++, iV += 2 )
			ublas::subrange( v, iV, iV + 2 ) = it->measurement;
	}


	/** creates a parameter vector from the network description */
	template< class VT >
	void buildParameterVector( VT& v )
	{
		namespace ublas = boost::numeric::ublas;
		std::size_t iV = 0;

		// free point positions
		for ( typename BundleAdjustmentNetwork< T >::PointsList::iterator it = m_net.points.begin(); 
			it != m_net.points.end(); it++, iV += 3 )
			ublas::subrange( v, iV, iV + 3 ) = *it;

		// body poses (the first defines the WCOS, and therefore is omitted)
		for ( typename BundleAdjustmentNetwork< T >::PoseList::iterator it = ++m_net.bodyPoses.begin();
			it != m_net.bodyPoses.end(); it++, iV += 6 )
		{
			ublas::vector_range< VT > trans( v, ublas::range( iV, iV + 3 ) );
			ublas::vector_range< VT > rot( v, ublas::range( iV + 3, iV + 6 ) );
			Math::vector_cast_assign( trans, it->translation() );
			Math::vector_cast_assign( rot, it->rotation().toLogarithm() );
			//TODO: it->rotation().toLogarithm( vr );
		}

		// camera extrinsics
		for ( typename BundleAdjustmentNetwork< T >::PoseList::iterator it = m_net.images.begin();
			it != m_net.images.end(); it++, iV += 6 )
		{
			ublas::vector_range< VT > trans( v, ublas::range( iV, iV + 3 ) );
			ublas::vector_range< VT > rot( v, ublas::range( iV + 3, iV + 6 ) );
			Math::vector_cast_assign( trans, it->translation() );
			Math::vector_cast_assign( rot, it->rotation().toLogarithm() );
			//TODO: it->rotation().toLogarithm( vr );
		}

		// camera intrinsics (pose + distortion)
		if ( m_net.bEstimateIntrinsics )
		{
			typename BundleAdjustmentNetwork< T >::IntrinsicList::iterator itI = m_net.intrinsics.begin();
			typename BundleAdjustmentNetwork< T >::DistortionList::iterator itD = m_net.distortions.begin();

			for ( ; itI != m_net.intrinsics.end(); itI++, itD++ )
			{
				v( iV++ ) = (*itI)( 0, 0 );
				v( iV++ ) = (*itI)( 0, 1 );
				v( iV++ ) = (*itI)( 0, 2 );
				v( iV++ ) = (*itI)( 1, 1 );
				v( iV++ ) = (*itI)( 1, 2 );

				ublas::subrange( v, iV, iV + 4 ) = *itD;
				iV += 4;
			}
		}
	}

	/** updates the parameters of the \c BundleAdjustmentNetwork from a given parameter vector */
	template< class VT >
	void updateParametersFromVector( VT& v )
	{
		namespace ublas = boost::numeric::ublas;
		std::size_t iV = 0;

		// free point positions
		for ( typename BundleAdjustmentNetwork< T >::PointsList::iterator it = m_net.points.begin(); 
			it != m_net.points.end(); it++, iV += 3 )
			*it = ublas::subrange( v, iV, iV + 3 );

		// body poses (the first defines the WCOS, and therefore is omitted)
		for ( typename BundleAdjustmentNetwork< T >::PoseList::iterator it = ++m_net.bodyPoses.begin();
			it != m_net.bodyPoses.end(); it++, iV += 6 )
		{
			*it = Math::Pose( 
				Math::Quaternion::fromLogarithm( ublas::subrange( v, iV + 3, iV + 6 ) ), 
				ublas::subrange( v, iV, iV + 3 ) );
		}

		// camera poses
		for ( typename BundleAdjustmentNetwork< T >::PoseList::iterator it = m_net.images.begin();
			it != m_net.images.end(); it++, iV += 6 )
		{
			*it = Math::Pose( 
				Math::Quaternion::fromLogarithm( ublas::subrange( v, iV + 3, iV + 6 ) ), 
				ublas::subrange( v, iV, iV + 3 ) );
		}

		// camera intrinsics (pose + distortion)
		if ( m_net.bEstimateIntrinsics )
		{
			typename BundleAdjustmentNetwork< T >::IntrinsicList::iterator itI = m_net.intrinsics.begin();
			typename BundleAdjustmentNetwork< T >::DistortionList::iterator itD = m_net.distortions.begin();

			for ( ; itI != m_net.intrinsics.end(); itI++, itD++, iV += 9 )
			{
				(*itI)( 0, 0 ) = v( iV + 0 );
				(*itI)( 0, 1 ) = v( iV + 1 );
				(*itI)( 0, 2 ) = v( iV + 2 );
				(*itI)( 1, 1 ) = v( iV + 3 );
				(*itI)( 1, 2 ) = v( iV + 4 );

				*itD = ublas::subrange( v, iV + 5, iV + 9 );
			}
		}
	}

protected:
	/**
	 * same as evaluateWithJacobian, but for a single 3D point in world coordinates
	 * @param result where to put the predicted 2d-measurement
	 * @param input the whole BA parameter vector
	 * @param J 2x<inputsize> jacobian of 2d-measurement wrt. all BA parameters
	 * @param pointJacobian 2x3 jacobian of 2D point wrt 3D input point
	 * @param camRotations precomputed camera rotation matrices (from input quaternions)
	 * @param bodyRotations precomputed body rotation matrices (from input quaternions)
	 * @param pIntrinsics precomputed or constant camera projection matrices
	 * @param iCamera index of camera intrinsics
	 * @param iImage index of camera pose
	 * @param p3d point in world coordinates
	 */
	template< class VT1, class VT2, class MT, class VT3 > 
	void evaluateSingleWorldPointWithJacobian( VT1& result, const VT2& input, MT& J, VT3& pointJacobian, 
		const std::vector< Math::Matrix< T, 3, 3 > >& camRotations,
		const std::vector< Math::Matrix< T, 3, 3 > >& bodyRotations,
		std::size_t iCamera, std::size_t iImage,
		const Math::Vector< T, 3 >& p3d ) const
	{
		namespace ublas = boost::numeric::ublas;
		// index into the parameter vector (used at various places in this routine)
		std::size_t iP;

		// transform point into camera coordinate frame
		iP = m_imageOffset + 6 * iImage;
		Math::Vector< T, 3 > transformed( ublas::prod( camRotations[ iImage ], p3d ) );
		noalias( transformed ) += ublas::subrange( input, iP, iP + 3 );

		// dehomogenize
		Math::Vector< T, 2 > dehomogenized;
		Function::Dehomogenization< 3 >().evaluate( dehomogenized, transformed );

		// apply intrinsics camera parameters
		Math::Vector< T, 2 > distorted;
		Math::Matrix< T, 2, 2 > jDistortion;
		if ( m_net.bEstimateIntrinsics )
		{
			// intrinsic parameters from parameter vector
			iP = m_intrinsicsOffset + 9 * iCamera;

			// distort
			Function::RadialDistortionWrtD< T >( dehomogenized )
				.evaluate( distorted, ublas::subrange( input, iP + 5, iP + 9 ) );

			// apply intrinsic matrix
			Function::CameraIntrinsicsMultiplication< T >( distorted )
				.evaluate( result, ublas::subrange( input, iP, iP + 5 ) );

			//std::cout << p3d << " -> " << result << std::endl;

			// jacobian wrt intrinsic matrix
			ublas::matrix_range< MT > jIntrinsics( J, ublas::range( 0, 2 ), ublas::range( iP, iP + 5 ) );
			Function::CameraIntrinsicsMultiplication< T >( distorted )
				.jacobian( ublas::subrange( input, iP, iP + 5 ), jIntrinsics );

			// jacobian of intrinsic multiplication wrt. p
			Math::Matrix< T, 2, 2 > K22;
			K22( 0, 0 ) = -input( iP ); K22( 0, 1 ) = -input( iP + 1 );
			K22( 1, 0 ) = 0;            K22( 1, 1 ) = -input( iP + 3 );

			// distortion jacobian wrt D
			Math::Matrix< T, 2, 4 > jDistD;
			Function::RadialDistortionWrtD< T >( dehomogenized )
				.jacobian( ublas::subrange( input, iP + 5, iP + 9 ), jDistD );
			noalias( ublas::subrange( J, 0, 2, iP + 5, iP + 9 ) ) = ublas::prod( K22, jDistD );

			// distortion jacobian wrt P
			Math::Vector< T, 4 > distCoeffs( ublas::subrange( input, iP + 5, iP + 5 + 4 ) );
			Math::Matrix< T, 2, 2 > jDistP;
			Function::RadialDistortionWrtP< T >( distCoeffs ).jacobian( dehomogenized, jDistP );
			noalias( jDistortion ) = ublas::prod( K22, jDistP );
		}
		else
		{
			// intrinsic parameters from network description

			// distort
			Function::RadialDistortionWrtD< T >( dehomogenized ).evaluate( distorted, m_net.distortions[ iCamera ] );

			// project
			const Math::Matrix< T, 3, 3 >& K( m_net.intrinsics[ iCamera ] );
			result( 0 ) = -( K( 0, 0 ) * distorted( 0 ) + K( 0, 1 ) * distorted( 1 ) + m_net.intrinsics[ iCamera ]( 0, 2 ) );
			result( 1 ) = -(                              K( 1, 1 ) * distorted( 1 ) + m_net.intrinsics[ iCamera ]( 1, 2 ) );

			// distortion jacobian wrt P
			Math::Matrix< T, 2, 2 > jDistP;
			Function::RadialDistortionWrtP< T >( m_net.distortions[ iCamera ] ).jacobian( dehomogenized, jDistP );
			noalias( jDistortion ) = -ublas::prod( ublas::subrange( K, 0, 2, 0, 2 ), jDistP );
		}

		// jacobian of dehomogenization
		Math::Matrix< T, 2, 3 > jDehom;
		Math::Matrix< T, 2, 3 > jDehom2;
		Function::Dehomogenization< 3 >().jacobian( transformed, jDehom2 );

		// multiply with dehomogenization jacobian
		noalias( jDehom ) = ublas::prod( jDistortion, jDehom2 );

		// jacobian wrt. camera translation
		iP = m_imageOffset + 6 * iImage;
		noalias( ublas::subrange( J, 0, 2, iP, iP + 3 ) ) = jDehom;

		// jacobian wrt. camera orientation
		Math::Matrix< T, 3, 3 > jCamOri;
		Function::LieRotation< T >( p3d ).jacobian( 
			ublas::subrange( input, iP + 3, iP + 6 ), jCamOri );
		noalias( ublas::subrange( J, 0, 2, iP + 3, iP + 6 ) ) = ublas::prod( jDehom, jCamOri );

		// jacobian
		noalias( pointJacobian ) = ublas::prod( jDehom, camRotations[ iImage ] );
	}

	// some offsets into the parameter vector
	
	/** offset of the first point in the parameter vector */
	std::size_t m_pointOffset;
	
	/** offset of the first body pose in the parameter vector */
	std::size_t m_bodyPoseOffset;
	
	/** offset of the first camera pose in the parameter vector */
	std::size_t m_imageOffset;
	
	/** offset of the first camera intrinsics value in the parameter vector */
	std::size_t m_intrinsicsOffset;
	

	/** reference to the network */
	BundleAdjustmentNetwork< T >& m_net;
};

} } // namespace Ubitrack::Calibration

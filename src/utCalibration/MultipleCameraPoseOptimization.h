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


#include <utMath/NewFunction/Function.h>
#include <utMath/NewFunction/Addition.h>
#include <utMath/NewFunction/Dehomogenization.h>
#include <utMath/NewFunction/LieRotation.h>
#include <utMath/NewFunction/LinearTransformation.h>
#include <utMeasurement/Measurement.h>

namespace Ubitrack { namespace Calibration {

#ifdef HAVE_LAPACK

/**
 * Function to minimize. Input is a 6-vector containing translation and exponential map rotation.
 */ 
template< class VType = double >
class ObjectiveFunction
{
public:
	ObjectiveFunction( const std::vector< Math::Vector< 3, VType > >& p3D, 
		const std::vector< Math::Matrix< 3, 3 > >& cameraRotations, 
		const std::vector< Math::Vector< 3 > >& cameraTranslations, 
		const std::vector< Math::Matrix< 3, 3, VType > >& cameraIntrinsics, 
		const std::vector< std::pair< unsigned, unsigned > > visibilities )
		: m_p3D( p3D )
		, m_camR( cameraRotations )
		, m_camT( cameraTranslations )
		, m_camI( cameraIntrinsics )
		, m_vis( visibilities )
	{}

	/**
	 * return the size of the result vector
	 */
	unsigned size() const
	{ return 2 * m_vis.size(); }


	/**
	 * @param result vector to store the result in
	 * @param input containing the parameters (target pose as 7-vector)
	 * @param J matrix to store the jacobian (evaluated for input) in
	 */
	template< class VT1, class VT2, class MT > 
	void evaluateWithJacobian( VT1& result, const VT2& input, MT& J ) const
	{
		namespace NF = Math::Function;
		namespace ublas = boost::numeric::ublas;
		
		for ( unsigned i = 0; i < m_vis.size(); i++ )
		{
				ublas::vector_range< VT1 > subResult( result, ublas::range( i * 2, ( i + 1 ) * 2 ) );
				ublas::matrix_range< MT > subJ( J, ublas::range( i * 2, ( i + 1 ) * 2 ), ublas::range( 0, 6 ) );
		
				( NF::Dehomogenization< 3 >() <<
				( NF::LinearTransformation< 3, 3 >( m_camI[ m_vis[ i ].second ] ) <<
					( NF::Addition< 3 >() <<
						( NF::fixedParameterRef< 3 >( m_camT[ m_vis[ i ].second ] ) ) <<
						( NF::LinearTransformation< 3, 3 >( m_camR[ m_vis[ i ].second ] ) <<
							( NF::Addition< 3 >() <<
								( NF::parameter< 3 >( 0 ) ) <<
								( NF::LieRotation() <<
									( NF::parameter< 3 >( 3  ) ) <<
									( NF::fixedParameterRef< 3 >( m_p3D[ m_vis[ i ].first ] ) )
								)
							)
						)
					)
				)
			).evaluateWithJacobian( input, subResult, subJ );
		}
	}
	
protected:
	const std::vector< Math::Vector< 3, VType > >& m_p3D;
	const std::vector< Math::Matrix< 3, 3 > >& m_camR;
	const std::vector< Math::Vector< 3 > >& m_camT;
	const std::vector< Math::Matrix< 3, 3, VType > >& m_camI;
	const std::vector< std::pair< unsigned, unsigned > > m_vis;
};


void checkConsistency (
	const std::vector < Math::Vector < 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector < 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar < double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< 3, 3 > >& camMatrices
	);

std::pair < Math::ErrorPose , double > 
	multipleCameraEstimatePose (
	const std::vector < Math::Vector < 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector < 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar < double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< 3, 3 > >& camMatrices,
	const int minCorrespondences,
	bool hasInitialPoseProvided,
	Math::Pose initialPose = Math::Pose(),
	int startIndex = 0,
	int endIndex = -1);

UBITRACK_EXPORT void multipleCameraPoseEstimationWithLocalBundles (
	const std::vector < Math::Vector < 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector < 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar < double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< 3, 3 > >& camMatrices,
	const int minCorrespondences,
	std::vector < Math::ErrorPose >& poses,
	std::vector < Math::Scalar < double > >& poseWeights,
	std::vector < Math::Scalar < int > >& localBundleSizes
	);

UBITRACK_EXPORT void multipleCameraPoseEstimation (
	const std::vector < Math::Vector < 3 > >&  points3d,
	const std::vector < std::vector < Math::Vector < 2 > > >& points2d,
	const std::vector < std::vector < Math::Scalar < double > > >& points2dWeights,
	const std::vector < Math::Pose >& camPoses,
	const std::vector < Math::Matrix< 3, 3 > >& camMatrices,
	const int minCorrespondences,
	Math::ErrorPose& pose,
	Math::Scalar < double > & poseWeight,
	bool hasInitialPoseProvided = false,
	Math::Pose initialPose = Math::Pose()
	);

#endif // HAVE_LAPACK

} } // namespace Ubitrack::Components

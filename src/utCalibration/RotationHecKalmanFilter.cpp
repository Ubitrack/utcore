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
 * Implements online rotation-only hand-eye-calibration
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#include "RotationHecKalmanFilter.h"
#ifdef HAVE_LAPACK
 
#include "Function/RotHecFunction.h"
#include <utTracking/Kalman.h>
#include <utMath/Function/VectorNormalize.h>
#include <utMath/CovarianceTransform.h>

namespace Ubitrack { namespace Calibration {

namespace ublas = boost::numeric::ublas;

RotationHecKalmanFilter::RotationHecKalmanFilter()
{
	m_state.value = Math::Vector< 4 >( 0, 0, 0, 1 );
	m_state.covariance = ublas::identity_matrix< double >( 4, 4 ) * 1e2;
}

void RotationHecKalmanFilter::addMeasurement( const Math::Quaternion& a, const Math::Quaternion& b )
{
	Math::ErrorVector< 4 > kalmanMeasurement;
	kalmanMeasurement.value = Math::Vector< 4, double >::zeros();
	kalmanMeasurement.covariance = ublas::identity_matrix< double >( 4 ) * 1e-2;
	// TODO: for error propagation use RotHecCombine

	// do the filter update
	Function::RotHecMeasurement mf( a, b.negateIfCloser( a ) );
	Tracking::kalmanMeasurementUpdate< 4, 4 >( m_state, mf, kalmanMeasurement, 0, m_state.value.size() );

	// normalize the result to ensure quaternion properties
	m_state = Math::transformWithCovariance< 4, 4 >( Math::Function::VectorNormalize( 4 ), m_state );
	
	m_state.covariance += ublas::identity_matrix< double >( 4 ) * 1e-12;
}

} }

#endif // HAVE_LAPACK

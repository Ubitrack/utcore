/**
 * @ingroup tracking_algorithms
 * @file
 * Functions for Correlation
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_CORRELATION_H_INCLUDED__
#define __UBITRACK_CALIBRATION_CORRELATION_H_INCLUDED__

#ifdef HAVE_LAPACK

#include <utCore.h>
#include <utMath/Pose.h>
#include <utMath/Vector.h>
#include <utMath/Scalar.h>
#include <vector>

#include <utUtil/Exception.h>

namespace Ubitrack { namespace Calibration {


UBITRACK_EXPORT double computeCorrelation ( const std::vector< Math::Vector3d >& left,
											const std::vector< Math::Vector3d >& right);


}};

#endif // HAVE_LAPACK

#endif

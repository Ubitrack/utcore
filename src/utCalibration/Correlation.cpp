
/**
 * @ingroup tracking_algorithms
 * @file
 * Implementation of Correlation
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */

#include "Correlation.h"
#ifdef HAVE_LAPACK

#include <utMath/Matrix.h>
#include <utUtil/Exception.h>

namespace Ubitrack { namespace Calibration {

double computeCorrelation ( const std::vector< Math::Vector< 3, double > >& left,
							const std::vector< Math::Vector< 3, double > >& right)
{
	return 1.0f;
}


} } // namespace Ubitrack::Calibration

#endif // HAVE_LAPACK

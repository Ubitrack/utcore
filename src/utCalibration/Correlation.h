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

#include <utCore.h> // UBITRACK_EXPORT
#include <vector> // std::vector

namespace Ubitrack { namespace Calibration {


UBITRACK_EXPORT double computeCorrelation ( const std::vector< double >& left,
											const std::vector< double >& right);


}};

#endif // HAVE_LAPACK

#endif

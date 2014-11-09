/**
 * @ingroup tracking_algorithms
 * @file
 * Functions for Correlation
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */

#ifndef __UBITRACK_ALGORITHM_CORRELATION_H_INCLUDED__
#define __UBITRACK_ALGORITHM_CORRELATION_H_INCLUDED__


#include <utCore.h> // UBITRACK_EXPORT

// std
#include <vector> // std::vector

namespace Ubitrack { namespace Algorithm {


UBITRACK_EXPORT double computeCorrelation ( const std::vector< double >& left,
											const std::vector< double >& right);

}}; // namespace Ubitrack::Algorithm

#endif // __UBITRACK_ALGORITHM_CORRELATION_H_INCLUDED__

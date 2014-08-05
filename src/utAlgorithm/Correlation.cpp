
/**
 * @ingroup tracking_algorithms
 * @file
 * Implementation of Correlation
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */

#include "Correlation.h"
#include <utMath/Stochastic/Correlation.h>
#include <math.h>



namespace Ubitrack { namespace Algorithm {


double computeCorrelation ( const std::vector< double >& left,
							const std::vector< double >& right)
{
	if( left.empty() && right.empty() )
		return 1;
	
	return Math::Stochastic::correlation( left.begin(), left.end(), right.begin(), right.end() );
}

} } // namespace Ubitrack::Algorithm


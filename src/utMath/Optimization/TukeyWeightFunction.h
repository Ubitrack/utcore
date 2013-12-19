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
 * @ingroup math
 * @file
 * Tukey weights for robust optimizers
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

# include <math.h>
#ifndef __UBITRACK_MATH_TUKEYWEIGHTFUNCTION_H_INCLUDED__
#define __UBITRACK_MATH_TUKEYWEIGHTFUNCTION_H_INCLUDED__

namespace Ubitrack { namespace Math { 

class TukeyWeightFunction
{
private:
	unsigned rowsPerMeasurement;
	double c;

public:
	TukeyWeightFunction( unsigned rowsPerMeasurement, double c )
	{
		this->c = c;
		this->rowsPerMeasurement = rowsPerMeasurement;
	}


	bool noWeights() const
	{ return false; }

	template< class VT1, class VT2 > 
	void computeWeights( const VT1& errorVector, VT2& weightVector ) const
	{
		for ( unsigned i = 0; i < errorVector.size(); i += rowsPerMeasurement )
		{	
			typename VT1::value_type ro, e = 0, weight;
			
			for ( unsigned j = 0; j < rowsPerMeasurement; j++ )
				e += errorVector( i + j ) * errorVector( i + j );

			if ( e != 0 )
			{
				if(!(e>c*c))
					ro=(pow(c,2)/6)*(1-pow((1- e/(c*c)),3));
				else
					ro=(c*c/6);
				weight = sqrt(ro)/ro;
			}
			else
				weight = 0;


			for ( unsigned j = 0; j < rowsPerMeasurement; j++ )
				weightVector( i + j ) = weight;
		}
	}
};

} } // namespace Ubitrack::Math

#endif

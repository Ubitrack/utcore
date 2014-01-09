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
 * function to read a pixel from an image
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */

#ifndef __UBITRACK_CALIBRATION_FUNCTION_IMAGELOOKUP_H_INCLUDED__
#define __UBITRACK_CALIBRATION_FUNCTION_IMAGELOOKUP_H_INCLUDED__

#include "MultiVariateFunction.h"
 
namespace Ubitrack { namespace Math { namespace Optimization { namespace Function {

/**
 * Function that reads a pixel from an image. The sobel operator is used to compute the derivative.
 */
template< class ImageT >
struct ImageLookup
	: public MultiVariateFunction< ImageLookup< ImageT >, 1 >
{
	/**
	 * Just stores the reference to the image
	 */
	ImageLookup( const Vision::Image& _image )
		: image( _image )
	{}
	
	/*
	 * @param result the pixel value
	 * @param p the pixel coordinate
	 */
	template< class VT1, class VT2 > 
	void evaluate( VT1& result, const VT2& p ) const
	{
		typedef typename VT1::value_type T;
		int x = static_cast< int >( std::floor( p( 0 ) ) );
		int y = static_cast< int >( std::floor( p( 1 ) ) );

		if ( x < 0 || x >= image.width - 1 || y < 0 || y >= image.height - 1 )
		{
			result( 0 ) = 0;
			return;
		}
		
		// do interpolation
		T dx = p( 0 ) - x;
		T dy = p( 1 ) - y;
		ImageT* i = reinterpret_cast< ImageT* >( image.imageData + y * image.widthStep ) + x;
		T a = i[ 0 ] + dx * ( i[ 1 ] - i[ 0 ] );
		i = reinterpret_cast< ImageT* >( reinterpret_cast< char* >( i ) + image.widthStep );
		T b = i[ 0 ] + dx * ( i[ 1 ] - i[ 0 ] );
		result( 0 ) = a + dy * ( b - a );
	}
	
	template< class LeftHand, class DestinationMatrix, class Param1 >
	void multiplyJacobian1( const LeftHand& l, DestinationMatrix& j, const Param1& p ) const
	{
		typedef typename LeftHand::value_type T;
		int x = static_cast< int >( std::floor( p( 0 ) ) );
		int y = static_cast< int >( std::floor( p( 1 ) ) );

		if ( p( 0 ) - T( x ) > T( 0.5 ) ) // round x
			x++;
		if ( p( 1 ) - T( y ) > T( 0.5 ) ) // round y
			y++;

		if ( x < 1 || x >= image.width - 2 || y < 1 || y >= image.height - 2 )
		{
			j( 0, 0 ) = T( 0 );
			j( 0, 1 ) = T( 0 );
			return;
		}
		
		// compute x and y sobel
		T sobX( 0 );
		T sobY( 0 );

		ImageT* i = reinterpret_cast< ImageT* >( image.imageData + ( y - 1 ) * image.widthStep ) + x - 1;
		sobX += T( i[ 2 ] ) - T( i[ 0 ] );
		sobY -= T( i[ 0 ] ) + 2 * T( i[ 1 ] ) + T( i[ 2 ] );
		
		i = reinterpret_cast< ImageT* >( reinterpret_cast< char* >( i ) + image.widthStep );
		sobX += 2 * ( T( i[ 2 ] ) - T( i[ 0 ] ) );
		
		i = reinterpret_cast< ImageT* >( reinterpret_cast< char* >( i ) + image.widthStep );
		sobX += T( i[ 2 ] ) - T( i[ 0 ] );
		sobY += T( i[ 0 ] ) + 2 * T( i[ 1 ] ) + T( i[ 2 ] );
		
		sobX /= 8;
		sobY /= 8;
		
		j( 0, 0 ) = l( 0, 0 ) * sobX;
		j( 0, 1 ) = l( 0, 0 ) * sobY;
	}
	
	/** reference to the image */
	const Vision::Image& image;
};

} } } // namespace Ubitrack::Calibration::Function

#endif

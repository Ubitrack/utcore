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
 * Class to store camera intrinsic parameters
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */


#ifndef __CAMERA_INTRINSICS_H_INCLUDED__
#define __CAMERA_INTRINSICS_H_INCLUDED__

//Boost
#include <boost/serialization/access.hpp>

//Ubitrack
#include <utMath/Matrix.h>
#include <utMath/Vector.h>
#include <utMath/Functors/MatrixFunctors.h>


namespace Ubitrack { namespace Math {

/**
 * @ingroup math
 * Stores all intrinsics camera parameters to have one compact 
 * representation for easier handling within dataflow configurations.
 * @param T can be double of float
 */

template < typename T > 
struct CameraIntrinsics
{
	/** definition of the matrix type */
	typedef Math::Matrix< 3, 3, T > matrix_type;
	
	/** definition of the radial distortion parameter type */
	typedef Math::Vector< 6, T > radial_type;
	
	/** definition of the tangential distortion parameter type */
	typedef Math::Vector< 2, T > tangential_type;
	
protected:
	/** Functor for the matrix inverse */
	const Math::Functors::matrix_inverse m_inverter;

public:		
	/**
	 * cameras' image calibration dimensions
	 * due to normalisation this should be 1, 1
	 */
	Ubitrack::Math::Vector< 2, std::size_t > dimension;
	
	/** cameras' 3x3-intrinsic matrix (normalized) */
	matrix_type matrix;
	
	/** inverse of cameras' 3x3-intrinsic matrix (normalized) */
	matrix_type matrix_inv;
	
	/** signs how many distortion parameters are actually used */
	std::size_t radial_size;
	
	/** radial distortion parameters */
	radial_type radial_params;
	
	/** tangential distortion parameters */
	tangential_type tangential_params;
		
	/** Standard Constructor */
	CameraIntrinsics(  )
		: dimension( Math::Vector< 2, std::size_t >( 1, 1 ) )
		, matrix( matrix_type::identity() )
		, matrix_inv( matrix_type::identity() )
		, radial_size( 0 )
		, radial_params( radial_type::zeros() )
		, tangential_params( tangential_type::zeros() )
		{}
	
	/** Constructor to use with old OpenCV values (2 radial distortion parameters) */	
	CameraIntrinsics( const Math::Matrix< 3, 3, T > &intrinsicMatrix, const Math::Vector< 2, T > &_radial, const Math::Vector< 2, T > &_tangential )
		: dimension( Math::Vector< 2, std::size_t >( 1, 1 ) )
		, matrix( intrinsicMatrix )
		, matrix_inv( m_inverter( intrinsicMatrix ) )
		, radial_size( 2 )
		, radial_params( radial_type::zeros() )
		, tangential_params( _tangential )
		{
			radial_params( 0 ) = _radial( 0 );
			radial_params( 1 ) = _radial( 1 );
		}
	
	/** Constructor to use with newer OpenCV values (6 radial distortion parameters) */	
	CameraIntrinsics( const Math::Matrix< 3, 3, T > &intrinsicMatrix, const Math::Vector< 6, T > &_radial, const Math::Vector< 2, T > &_tangential )
		: dimension( Math::Vector< 2, std::size_t >( 1, 1 ) )
		, matrix( intrinsicMatrix )
		, matrix_inv( m_inverter( intrinsicMatrix ) )
		, radial_size( 6 )
		, radial_params( _radial )
		, tangential_params( _tangential )
		{
		}

protected:
	friend class ::boost::serialization::access;
	
	/** boost::serialization helper */
	template< class Archive >
	void serialize( Archive& ar, const unsigned int version )
	{
		//the casts are necessary to set const values...
		for( std::size_t i = 0; i < 9; ++i )
			ar & matrix( ( i / 3 ), ( i % 3 ) );
		
		// only necessary when loading:
		matrix_inv = m_inverter( matrix );
		
		ar & dimension( 0 );
		ar & dimension( 1 );		
		ar & radial_size );
		
		for( std::size_t i = 0; i < radial_size; ++i )
			ar & radial_params[ i ];
		
		std::size_t tan_dim = 2;
		ar & tan_dim;
		
		ar & tangential_params[ 0 ];
		ar & tangential_params[ 1 ];
	}
};

} } // namespace Ubitrack::Math

#endif //__CAMERA_INTRINSICS_H_INCLUDED__


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


#ifndef __UBITRACK_MATH_CAMERA_INTRINSICS_H_INCLUDED__
#define __UBITRACK_MATH_CAMERA_INTRINSICS_H_INCLUDED__

//std
#include <iosfwd>	// std::ostream
#include <math.h>	// std::atan, std::acos

//Boost
#include <boost/serialization/access.hpp>

//Ubitrack
#include <utMath/Vector.h>
#include <utMath/Matrix.h>

namespace Ubitrack { namespace Math {

/**
 * @ingroup math
 * Stores all intrinsics camera parameters to have one compact 
 * representation for easier handling within dataflow configurations
 * and algorithms.
 * @param T can be \c double of \c float
 */

template < typename T > 
struct CameraIntrinsics
{
	/** defining the own type */
	typedef CameraIntrinsics		self_type;
	
	/** definition of the matrix type */
	typedef Math::Matrix< T, 3, 3 >	matrix_type;
	
	/** definition of the radial distortion parameter type */
	typedef Math::Vector< T, 6 >	radial_type;
	
	/** definition of the tangential distortion parameter type */
	typedef Math::Vector< T, 2 >	tangential_type;
		
	/// describes all various calibration types offered so far:
	/// IMPLEMENTATION_RADIAL_TANGENTIAL[_SPECIALIZATION]
	enum CalibType { UNKNOWN, OPENCV_2_2, OPENCV_3_2, OPENCV_6_2, OPENCV_4_0_FISHEYE };
	
	/// provides the number of elements for the radial distortion ( use with \c CalibType )
	static const std::size_t radial( const int idx )
	{
		static const std::size_t a[] = { 0u, 2u, 3u, 6u, 4u };
		return a[ idx ];
	}
	
	/// provides the number of elements for the tangential distortion ( use with \c CalibType )
	static const std::size_t tangential( const int idx )
	{
		static const std::size_t a[] = { 0u, 2u, 2u, 2u, 0u };
		return a[ idx ];
	}

public:	

	/** defining the type/version of the calibration internally. */
	CalibType calib_type;
	
	/**
	 * cameras' image calibration dimensions
	 * if it is 0 then the dimension is unknown
	 */
	Ubitrack::Math::Vector< std::size_t, 2 > dimension;
	
	/** cameras' 3x3-intrinsic matrix (normalized) */
	matrix_type matrix;
	
	/// maybe switch to this kind of structure
	// struct Distortion
	// {
			
		// std::size_t size_radial;
		// Math::Vector< T, 6 > radial;
		
		// std::size_t size_tangetial;
		// Math::Vector< T, 2 > tangential;
		// Distortion( const Math::Vector< T, 6 >& rad, const Math::Vector< T, 2 > tan )
			// : raidal( rad )
			// , tangential( tan )
			// {}
	// };
	
	// Distortion distortion;
	
	/** signs how many distortion parameters are actually used */
	std::size_t radial_size;
	
	/** radial distortion parameters */
	radial_type radial_params;
	
	/** tangential distortion parameters */
	tangential_type tangential_params;
		
	/** Standard Constructor */
	CameraIntrinsics(  )
		: calib_type( UNKNOWN )
		, dimension( Math::Vector< std::size_t, 2 >( 0, 0 ) )
		, matrix( matrix_type::identity() )
		, radial_size( 0 )
		, radial_params( radial_type::zeros() )
		, tangential_params( tangential_type::zeros() )
		{}
	
	/** Constructor to use with old OpenCV values (2 radial distortion parameters) */	
	CameraIntrinsics(const Math::Matrix< T, 3, 3 > &intrinsicMatrix, const Math::Vector< T, 2 > &_radial, const Math::Vector< T, 2 > &_tangential, const std::size_t width = 0, const std::size_t height = 0)
		: calib_type( OPENCV_2_2 )
		, dimension( Math::Vector< std::size_t, 2 >( width, height ) )
		, matrix( intrinsicMatrix )
		, radial_size( 2 )
		, radial_params( radial_type::zeros() )
		, tangential_params( _tangential )
		{
			radial_params( 0 ) = _radial( 0 );
			radial_params( 1 ) = _radial( 1 );
		}
		
	/** Constructor to use with newer OpenCV values (3 radial distortion parameters) */	
	CameraIntrinsics(const Math::Matrix< T, 3, 3 > &intrinsicMatrix, const Math::Vector< T, 3 > &_radial, const Math::Vector< T, 2 > &_tangential, const std::size_t width = 0, const std::size_t height = 0)
		: calib_type( OPENCV_3_2 )
		, dimension(Math::Vector< std::size_t, 2 >(width, height))
		, matrix( intrinsicMatrix )
		, radial_size( 2 )
		, radial_params( radial_type::zeros() )
		, tangential_params( _tangential )
		{
			radial_params( 0 ) = _radial( 0 );
			radial_params( 1 ) = _radial( 1 );
			radial_params( 2 ) = _radial( 2 );
		}
	
	/** Constructor to use with newer OpenCV values (6 radial distortion parameters) */	
	CameraIntrinsics(const Math::Matrix< T, 3, 3 > &intrinsicMatrix, const Math::Vector< T, 6 > &_radial, const Math::Vector< T, 2 > &_tangential, const std::size_t width = 0, const std::size_t height = 0)
		: calib_type( OPENCV_6_2 )
		, dimension(Math::Vector< std::size_t, 2 >(width, height))
		, matrix( intrinsicMatrix )
		, radial_size( 6 )
		, radial_params( _radial )
		, tangential_params( _tangential )
		{
		}
		
	/** Constructor to use with new OpenCV fish-eye values (4 distortion parameters, not tangential) */	
	CameraIntrinsics( const Math::Matrix< T, 3, 3 > &intrinsicMatrix, const Math::Vector< T, 4 > &_radial )
		: calib_type( OPENCV_4_0_FISHEYE )
		, dimension(Math::Vector< std::size_t, 2 >(width, height))
		, matrix( intrinsicMatrix )
		, radial_size( 4 )
		, tangential_params( tangential_type::zeros() )
		{
			radial_params( 0 ) = _radial( 0 );
			radial_params( 1 ) = _radial( 1 );
			radial_params( 2 ) = _radial( 2 );
			radial_params( 3 ) = _radial( 3 );
		}
	
	/** Constructor for the very general use case, accepting everything important */
	template< typename IntType >
	CameraIntrinsics(const Math::Vector< IntType, 2 >& size, const Math::Matrix< T, 3, 3 > &intrinsicMatrix, const std::size_t radSize, const Math::Vector< T, 6 > &_radial, const Math::Vector< T, 2 > &_tangential, const std::size_t width = 0, const std::size_t height = 0)
		: calib_type( UNKNOWN )
		, dimension(Math::Vector< std::size_t, 2 >(width, height))
		, matrix( intrinsicMatrix )
		, radial_size( radSize )
		, radial_params( _radial )
		, tangential_params( _tangential )
		{
			reset();
		}
		
	CameraIntrinsics& operator= ( const CameraIntrinsics< T >& rhs )
	{
		///@todo maybe use a swap function
		this->calib_type		= rhs.calib_type;
		this->dimension			= rhs.dimension;
		this->matrix			= rhs.matrix;
		this->radial_size		= rhs.radial_size;
		this->radial_params		= rhs.radial_params;
		this->tangential_params	= rhs.tangential_params;
		return *this;
	}
	
	void reset()
	{
		switch( radial_size )
		{
			case 2 : // old school calibration
				calib_type = OPENCV_2_2;
				break;
			case 3 : // new standard calibration
				calib_type = OPENCV_3_2;
				break;
			case 4: // fish-eye parameters
				calib_type = OPENCV_4_0_FISHEYE;
				break;
			case 6 : // wide-angle parameters
				calib_type = OPENCV_6_2;
				break;
		}
	}
	
	/** returns the opening angle of the lens in vertical direction */
	T angleVertical( ) const
	{
		// angle ~= 2*arctan( d/(2*f)) -> rad -> degree: * (180/pi) 
		static const T pi = std::acos( static_cast< T >( -1 ) ); //pi ~= acos(-1)
		static const T factor = (360/pi);
		return factor * std::atan( dimension( 1 ) / ( 2* matrix( 1, 1 ) ) );
	}
	
	/** returns the opening angle of the lens in horizontal direction */
	T angleHorizontal( ) const
	{
		// angle ~= 2*arctan( d/(2*f)) -> rad -> degree: * (180/pi) 
		static const T pi = std::acos( static_cast< T >( -1 ) ); //pi ~= acos(-1)
		static const T factor = (360/pi);
		return factor * std::atan( dimension( 0 ) / ( 2* matrix( 0, 0 ) ) );
	}
	
	/** returns the opening angle of the lens in diagonal direction */
	T angleDiagonal( ) const
	{
		// angle ~= 2*arctan( d/(2*f)) -> rad -> degree: * (180/pi) 
		static const T pi = std::acos( static_cast< T >( -1 ) ); //pi ~= acos(-1)
		static const T factor = (360/pi);
		const std::size_t squaredDiag  = ( dimension( 0 )*dimension( 0 )+dimension( 1 )*dimension( 1 ) );
		const T squaredFocal = ( matrix( 0, 0 )*matrix( 0, 0 )+matrix( 1, 1 )*matrix( 1, 1 ) );
		return factor * std::atan( std::sqrt( squaredDiag / (4*squaredFocal) ) );
	}
	
	/// corrects the Ubitrack intrinsics matrix to the corresponding left-handed (e.g. \c OpenCV ) camera matrix
	template< typename PrecisionType >	
	inline void flipHandiness( Math::CameraIntrinsics< PrecisionType >& intrinsics )
	{
		// compensate for left-handed OpenCV coordinate frame
		intrinsics.matrix ( 0, 2 ) *= -1;
		intrinsics.matrix ( 1, 2 ) *= -1;
		intrinsics.matrix ( 2, 2 ) *= -1;
	}
	

protected:
	friend class ::boost::serialization::access;
	
	/** boost::serialization helper */
	template< class Archive >
	void serialize( Archive& ar, const unsigned int version )
	{
		for( std::size_t i = 0; i < 9; ++i )
			ar & matrix( ( i / 3 ), ( i % 3 ) );
		
		ar & dimension( 0 );
		ar & dimension( 1 );		
		ar & radial_size;
		
		for( std::size_t i = 0; i < radial_size; ++i )
			ar & radial_params[ i ];
		
		std::size_t tan_dim = 2;
		ar & tan_dim;
		
		ar & tangential_params[ 0 ];
		ar & tangential_params[ 1 ];
		
		reset();
	}
};


/** @internal stream output operator */
template< typename T >
std::ostream& operator<<( std::ostream& s, const CameraIntrinsics< T >& intrCam )
{
	s << "Intrinsic matrix:\n";
	s << intrCam.matrix;
	s << "Resolution [width x height] : [" << intrCam.dimension( 0 ) << " x " << intrCam.dimension( 1 )<< "]";
	s << "\nAppr. field-of-view [HxVxD] : [" <<  intrCam.angleHorizontal() << " x " << intrCam.angleVertical() << " x " << intrCam.angleDiagonal() << "]";
	s << "\nTangential distortion   (2) : [" << intrCam.tangential_params( 0 ) << ", " << intrCam.tangential_params( 1 ) << "]";
	if( intrCam.radial_size )
	{
		s << "\nRadial distortion       (" << intrCam.radial_size << ") : [" << intrCam.radial_params( 0 );
		for( std::size_t i( 1 ); i<intrCam.radial_size; ++i )
			s <<  ", " << intrCam.radial_params[ i ];
		s << "]";
	}
	return s;
}

} } // namespace Ubitrack::Math

#endif //__UBITRACK_MATH_CAMERA_INTRINSICS_H_INCLUDED__


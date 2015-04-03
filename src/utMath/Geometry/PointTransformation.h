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
 * @ingroup math geometry
 * @file
 * Functions to spatially transform a 3D point.
 *
 * @author Christian Waechter <christian.waechter@in.tum.de>
 */ 

 
 
#ifndef __H__POINT_TRANSFORMATION_FUNCTIONS__
#define __H__POINT_TRANSFORMATION_FUNCTIONS__


#include <utMath/Vector.h> //includes static assert
#include <utMath/Matrix.h>


#include "container_traits.h"
#include "../Util/type_traits.h"
#include "../Stochastic/identity_iterator.h"

#include <algorithm> //std::transform

namespace Ubitrack { namespace Math { namespace Geometry {


/**
 * @ingroup math geometry functor
 * @brief Functor to spatially transform a \b 2D or \b 3D \b point .
 *
 * This functors provides several overloaded bracket operators (\c operator() )for
 * different types of \b 3D \b point and transformation matrix representation.
 *
 * Possible transformation matrix representations:
 * - \b Matrix2x3 (for \b 2D \b point transformation)
 * - \b Matrix3x3 (for \b 2D \b point transformation(e.g. homogeneous representation) or \b 3D \b point transformation )
 * - \b Matrix3x4 (for \b 3D \b point transformation)
 * - \b Matrix4x4 (for \b 3D \b point transformation, homogeneous representation)
 *
 * Possible \b 2D \b point representations:
 * - \b Vector2 (3rd dimension is assumed as 1)
 * - \b Vector3 (common representation)
 *
 * Possible \b 3D \b point representations:
 * - \b Vector2 (homogeneous point, 3rd dimension is assumed as 0)
 * - \b Vector3 (common representation)
 * - \b Vector4 (e.g. homogeneous representation with one for last dimension, but can also be different )
 */
struct TransformPoint
{
public:
	// internal template to catch wrong vector types and print an error message
	template< typename notSupportedMatrixType, typename notSupportedVectorType >
	notSupportedVectorType operator() ( const notSupportedMatrixType &, const notSupportedVectorType &  ) const
	{
		UBITRACK_STATIC_ASSERT( false, USE_ONLY_SUPPORTED_VECTOR_MATRIX_TYPES );
		return notSupportedVectorType();
	}
	
	/// @internal Specialization of bracket operator (\c operator() ) for \b 2-by-3 \b transformation of \b 2D \b points (as Vector2D)
	template< typename T >
	Math::Vector< T, 2 > operator() ( const Math::Matrix< T, 2, 3 > &transMat, const Math::Vector< T, 2 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 );
		return Math::Vector< T, 2 > ( e1, e2 );
	}
	
	/// @internal Specialization of bracket operator (\c operator() )for \b 2-by-3 \b transformation of \b 2D \b points (as Vector3D)
	template< typename T >
	Math::Vector< T, 2 > operator() ( const Math::Matrix< T, 2, 3 > &transMat, const Math::Vector< T, 3 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 ) * vec( 2 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 ) * vec( 2 );
		return Math::Vector< T, 2 > ( e1, e2 );
	}

	/// @internal Specialization of bracket operator (\c operator() ) for \b 3-by-3 \b transformation of \b 2D \b points (as Vector2D)
	template< typename T >
	Math::Vector< T, 3 > operator() ( const Math::Matrix< T, 3, 3 > &transMat, const Math::Vector< T, 2 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 2 );
		return Math::Vector< T, 3 > ( e1, e2, e3 );
	}
	
	/// @internal Specialization of bracket operator (\c operator() ) for \b 3-by-3 \b transformation of \b 2D \b points (as Vector3D)
	template< typename T >
	Math::Vector< T, 3 > operator() ( const Math::Matrix< T, 3, 3 > &transMat, const Math::Vector< T, 3 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 ) * vec( 2 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 ) * vec( 2 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 2 ) * vec( 2 );
		return Math::Vector< T, 3 > ( e1, e2, e3 );
	}

	/// @internal Specialization of bracket operator (\c operator() ) for \b 3-by-4 \b transformation of \b 3D \b points (as Vector2D)
	template< typename T >
	Math::Vector< T, 3 > operator() ( const Math::Matrix< T, 3, 4 > &transMat, const Math::Vector< T, 2 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 3 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 3 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 3 );
		return Math::Vector< T, 3 > ( e1, e2, e3 );
	}
	
	/// @internal Specialization of bracket operator (\c operator() ) for \b 3-by-4 \b transformation of \b 3D \b points (as Vector3D)
	template< typename T >
	Math::Vector< T, 3 > operator() ( const Math::Matrix< T, 3, 4 > &transMat, const Math::Vector< T, 3 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 ) * vec( 2 ) + transMat( 0, 3 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 ) * vec( 2 ) + transMat( 1, 3 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 2 ) * vec( 2 ) + transMat( 2, 3 );
		return Math::Vector< T, 3 > ( e1, e2, e3 );
	}
	
	/// @internal Specialization of bracket operator (\c operator() ) for \b 3-by-4 \b transformation of \b 3D \b points (as 4D vector)
	template< typename T >
	Math::Vector< T, 3 > operator() ( const Math::Matrix< T, 3, 4 > &transMat, const Math::Vector< T, 4 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 ) * vec( 2 ) + transMat( 0, 3 ) * vec( 3 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 ) * vec( 2 ) + transMat( 1, 3 ) * vec( 3 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 2 ) * vec( 2 ) + transMat( 2, 3 ) * vec( 3 );
		return Math::Vector< T, 3 > ( e1, e2, e3 );
	}
	
	/// @internal Specialization of bracket operator (\c operator() ) for \b 4-by-4 \b transformation of \b 2D \b points
	template< typename T >
	Math::Vector< T, 4 > operator() ( const Math::Matrix< T, 4, 4 > &transMat, const Math::Vector< T, 2 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 3 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 3 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 3 );
		const T e4 = transMat( 3, 0 ) * vec( 0 ) + transMat( 3, 1 ) * vec( 1 ) + transMat( 3, 3 );
		return Math::Vector< T, 4 > ( e1, e2, e3, e4 );
	}
	
	/// @internal Specialization of bracket operator (\c operator() ) for \b 4-by-4 \b transformation of \b 3D \b points
	template< typename T >
	Math::Vector< T, 4 > operator() ( const Math::Matrix< T, 4, 4 > &transMat, const Math::Vector< T, 3 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 ) * vec( 2 ) + transMat( 0, 3 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 ) * vec( 2 ) + transMat( 1, 3 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 2 ) * vec( 2 ) + transMat( 2, 3 );
		const T e4 = transMat( 3, 0 ) * vec( 0 ) + transMat( 3, 1 ) * vec( 1 ) + transMat( 3, 2 ) * vec( 2 ) + transMat( 3, 3 );
		return Math::Vector< T, 4 > ( e1, e2, e3, e4 );
	}
	
	/// @internal Specialization of bracket operator (\c operator() ) for \b 4-by-4 \b transformation of \b 4D \b points
	template< typename T >
	Math::Vector< T, 4 > operator() ( const Math::Matrix< T, 4, 4 > &transMat, const Math::Vector< T, 4 > &vec ) const
	{
		const T e1 = transMat( 0, 0 ) * vec( 0 ) + transMat( 0, 1 ) * vec( 1 ) + transMat( 0, 2 ) * vec( 2 ) + transMat( 0, 3 ) * vec( 3 );
		const T e2 = transMat( 1, 0 ) * vec( 0 ) + transMat( 1, 1 ) * vec( 1 ) + transMat( 1, 2 ) * vec( 2 ) + transMat( 1, 3 ) * vec( 3 );
		const T e3 = transMat( 2, 0 ) * vec( 0 ) + transMat( 2, 1 ) * vec( 1 ) + transMat( 2, 2 ) * vec( 2 ) + transMat( 2, 3 ) * vec( 3 );
		const T e4 = transMat( 3, 0 ) * vec( 0 ) + transMat( 3, 1 ) * vec( 1 ) + transMat( 3, 3 ) * vec( 2 ) + transMat( 3, 3 ) * vec( 3 );
		return Math::Vector< T, 4 > ( e1, e2, e3, e4 );
	}
};


/**
 * @ingroup math geometry
 * @brief transforms several points spatially using iterators pointing to the storage class of the points.
 *
 * This function can be applied to nearly any storage class providing 
 * access to the single points via forward iterators
 * (e.g. \c std::vector, \c std::list or \c std::set, etc ).
 * Although designed for stl-containers it is not limited to these ones
 * and this function can be applied to fixed size arrays as well or other storage classes.
 *
 * The function can transform either \b 2D, \b 3D or \b 4D points and therefore assumes homogeneous 
 * coordinates for the lower dimensional cases ( \b 2D and \b 3D ).
 * It can perform the following actions:
 * - \b 2D : @f$ \hat{p}_{2x1} = M_{2x3} \cdot [p_{1} p_{2} 1]^T @f$ 
 * - \b 2D : @f$ \hat{p}_{2x1} = M_{2x3} \cdot [p_{1} p_{2} p_{3}]^T @f$ 
 * - \b 2D : @f$ \hat{p}_{3x1} = M_{3x3} \cdot [p_{1} p_{2} 1]^T @f$ 
 * - \b 2D : @f$ \hat{p}_{3x1} = M_{3x3} \cdot [p_{1} p_{2} p_{3}]^T @f$
 * - \b 3D : @f$ \hat{p}_{3x1} = M_{3x4} \cdot [p_{1} p_{2} 0 1]^T @f$ 
 * - \b 3D : @f$ \hat{p}_{3x1} = M_{3x4} \cdot [p_{1} p_{2} p_{3} 1]^T @f$ 
 * - \b 3D : @f$ \hat{p}_{3x1} = M_{3x4} \cdot [p_{1} p_{2} p_{3} p_{4}]^T @f$ 
 * - \b 3D : @f$ \hat{p}_{4x1} = M_{4x4} \cdot [p_{1} p_{2} 0 1]^T @f$
 * - \b 3D : @f$ \hat{p}_{4x1} = M_{4x4} \cdot [p_{1} p_{2} p_{3} 1]^T @f$
 * - \b 3D : @f$ \hat{p}_{4x1} = M_{4x4} \cdot [p_{1} p_{2} p_{3} p_{4}]^T @f$
 * 
 * Example use case:\n
 @code
 Matrix< double, 3, 4 > trans; // <- should be filled with values
 std::vector< Vector3d > points3d; // <- should be filled with values
 transform_points( trans, points3d.begin(), points3d.end(), points3d.begin() );
 // or
 std::vector< Vector3d > points3dOut; // <- will be filled with values
 points3dOut.reserve( n ); // <- storage allocation via reserve() and number of elements( =n )
 transform_points( trans, points3d.begin(), points3d.end(), std::back_inserter( points3dOut ) ); 
 @endcode
 * 
 * @tparam T built-in type of matrix and input/output vectors ( e.g. \c double or \c float )
 * @tparam M first dimension of matrix (rows)
 * @tparam N second dimension of matrix (columns)
 * @tparam ForwardIterator1 type of forward iterator to container of input points
 * @tparam ForwardIterator2 type of the output iterator for the projected points
 * @param transformation the \b 3-by-4 \b transformation \b matrix applied to spatially transform the \b points
 * @param iBegin \c iterator pointing to first element in the input container/storage class of the \b points ( usually \c begin() )
 * @param iEnd \c iterator pointing behind the last element in the input container/storage class of the \b points ( usually \c end() )
 * @param iOut output \c iterator pointing to first element in container/storage class for storing the transformed points as a result ( usually \c begin() or \c std::back_inserter(container) )
 */
template< typename T, std::size_t M, std::size_t N, typename ForwardIterator1, typename ForwardIterator2 >
inline void transform_points( const Math::Matrix< T, M, N > &transformation, const ForwardIterator1 iBegin, const ForwardIterator1 iEnd, ForwardIterator2 iOut )
{

	// determine the types of the iterators.
	// Since output iterators (e.g. std::back_insert_itertator via std::back_inserter)
	// can be used as well a simple value_type is not enough
	// -> therefore designed an own container_traits struct
	typedef typename Ubitrack::Util::container_traits< ForwardIterator1 >::value_type vector_type_in;
	typedef typename Ubitrack::Util::container_traits< ForwardIterator2 >::value_type vector_type_out;
	typedef typename vector_type_in::value_type value_type_in;
	typedef typename vector_type_out::value_type value_type_out;
		
	UBITRACK_STATIC_ASSERT( (((M == 2) && (N == 3)) || ((M == 3) && (N == 3)) || ((M == 3) && (N == 4)) || ((M == 4) && (N == 4))), USING_A_NON_STANDARD_TRANSFORMATION_MATRIX );
	UBITRACK_STATIC_ASSERT( (Ubitrack::Util::is_same< value_type_in, T >::value ), MATRIX_AND_VECTORS_NEED_SAME_BUILTIN_TYPE ); // e.g. only float or only double
	UBITRACK_STATIC_ASSERT( (Ubitrack::Util::is_same< value_type_in, value_type_out >::value ), INPUT_AND_OUTPUT_VECTOR_NEED_SAME_BUILTIN_TYPE );
	UBITRACK_STATIC_ASSERT( (Ubitrack::Util::is_same< vector_type_out, Math::Vector< T, M > >::value ), OUTPUT_VECTOR_NEEDS_SAME_DIMENSION_AS_MATRIX_ROWS );

	// std::transform( iBegin, iEnd, iOut, std::bind1st( TransformPoint< T, M, N, vector_type_in >(), transformation ) );
	const std::size_t n = std::distance( iBegin, iEnd );
	Ubitrack::Util::identity< const Math::Matrix< T, M, N > > id_container( transformation, n );
	std::transform( id_container.begin(), id_container.end(), iBegin, iOut, TransformPoint() );
}


} } } // namespace Ubitrack::Math::Geometry

#endif //__H__POINT_TRANSFORMATION_FUNCTIONS__



#ifndef __UBITRACK_MATH_UTIL_MATRIX_TRAITS_H_INCLUDED__
#define __UBITRACK_MATH_UTIL_MATRIX_TRAITS_H_INCLUDED__

#include "type_traits.h"
//#include <type_traits> //integral_constant <- commented to compile on standard gcc (without c++0x features)

#include <cstddef> // std::size_t

namespace Ubitrack { namespace Math {

//some forward declaration
template< typename T, std::size_t M, std::size_t N > class Matrix;

namespace Util {


struct fixed_matrix_storage_tag {};
struct dynamic_matrix_storage_tag {};
struct unknown_matrix_storage_tag {};

/** @internal specialization for unbounded Matrices */
template< typename T >
struct has_fixed_storage< Math::Matrix< T, 0, 0 > >
	: public Ubitrack::Util::false_type{};

/** @internal specialization for bounded Matrices */
template< typename T, std::size_t M, std::size_t N >
struct has_fixed_storage< Math::Matrix< T, M, N > >
	: public Ubitrack::Util::true_type{};

/** @internal specialization for unbounded Matrices */
template< typename T >
struct has_dynamic_storage< Math::Matrix< T, 0, 0 > >
	: public Ubitrack::Util::true_type{};

/** @internal specialization for bounded Matrices */
template< typename T, std::size_t M, std::size_t N >
struct has_dynamic_storage< Math::Matrix< T, M, N > >
	: public Ubitrack::Util::false_type{};



/**
 * @internal MatrixTraits
 *
 * This functor supports the writing of compile time dependend code
 * by providing general type information of the underlaying
 * Matrix. Due to templating this functor can be overloaded for
 * various types of matrix by template specialization.
 */
template< typename MatType >
struct matrix_traits
{
	typedef MatType self_type;
	typedef unknown_storage_tag storage_category;
	typedef unknown_matrix_storage_tag storage_type;
	typedef typename std::size_t size_type;
	typedef typename MatType::value_type value_type;
};

/** @internal specialization for unbounded Matrices */
template< typename T >
struct matrix_traits< Math::Matrix< T, 0, 0 > >
{
	// typedef typename Math::Vector< T, 0 >::base_type base_type;
	typedef typename Math::Matrix< T, 0, 0 > self_type;
	typedef dynamic_storage_tag storage_category;
	typedef dynamic_matrix_storage_tag storage_type;
	
	typedef typename std::size_t size_type;
	typedef T value_type;
	
	// names a function to determine the size of a vector dynamically
	static const size_type size1( const self_type& rhs ) 
	{
		return rhs.size1();
	}
	
	static const size_type size2( const self_type& rhs ) 
	{
		return rhs.size2();
	}
	
	static T* ptr( const self_type& rhs )
	{
		return const_cast< T* > ( rhs.data().begin() );
	}
};

/** @internal specialization for bounded Matrices */
template< typename T, std::size_t M, std::size_t N >
struct matrix_traits< Math::Matrix< T, M, N > >
{
	// typedef typename Math::Vector< T, N >::base_type base_type;
	typedef typename Math::Matrix< T, M, N > self_type;
	typedef fixed_storage_tag storage_category;
	typedef fixed_matrix_storage_tag storage_type;
	typedef typename std::size_t size_type;
	typedef T value_type;

	static const size_type size1 = M;
	static const size_type size2 = N;
	
	static T* ptr( const self_type& rhs )
	{
		return const_cast< T* > ( rhs.data().begin() );
	}
	
};

} } } // namespace Ubitrack::Math::Util

#endif //__UBITRACK_MATH_UTIL_MATRIX_TRAITS_H_INCLUDED__

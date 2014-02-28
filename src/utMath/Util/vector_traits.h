

#ifndef __UBITRACK_MATH_UTIL_VECTOR_TRAITS_H_INCLUDED__
#define __UBITRACK_MATH_UTIL_VECTOR_TRAITS_H_INCLUDED__


//#include <type_traits> //integral_constant <- commented to compile on standard gcc (without c++0x features)
#include <cstddef> // std::size_t

namespace Ubitrack { namespace Math {

//some forward declaration
template< typename T, std::size_t N > class Vector;

namespace Util {

/**
 * @internal Constant value struct
 *
 * This functor reassembles the integral constant struct,
 * similar to the c++0x standard and supports the writing
 * of compile time dependend code.
 * see also http://en.cppreference.com/w/cpp/types/integral_constant
 */
template< class T, T v >
struct constant_value
{
	static const T value = v;
	typedef T value_type;
	typedef constant_value< T, v > type;
	// operator T() { return v; } //<- automatic type conversion, not activated right now
};

/** @internal specialization for boolean true type */
typedef constant_value< bool, true > true_type;

/** @internal specialization for boolean false type */
typedef constant_value< bool, false > false_type;

/** @internal identifies if a Vectors' size is known at compile time */
template< typename VType >
struct has_fixed_storage
	: public false_type{};

/** @internal specialization for unbounded Vectors */
template< typename T >
struct has_fixed_storage< Math::Vector< T, 0 > >
	: public false_type{};

/** @internal specialization for bounded Vectors */
template< typename T, std::size_t N >
struct has_fixed_storage< Math::Vector< T, N > >
	: public true_type{};

/** @internal identifies if a Vectors' size is unknown at compile time */
template< typename VType >
struct has_dynamic_storage
	: public true_type{};

/** @internal specialization for unbounded Vectors */
template< typename T >
struct has_dynamic_storage< Math::Vector< T, 0 > >
	: public true_type{};

/** @internal specialization for bounded Vectors */
template< typename T, std::size_t N >
struct has_dynamic_storage< Math::Vector< T, N > >
	: public false_type{};


struct fixed_storage_tag {};
struct dynamic_storage_tag {};
struct unknown_storage_tag {};


/**
 * @internal Vector_Traits
 *
 * This functor supports the writing of compile time dependend code
 * by providing general type information of the underlaying
 * Vector. Due to templating this functor can be overloaded for
 * various types of vector by template specialization.
 */
template< typename VType >
struct vector_traits
{
	typedef VType self_type;
	typedef unknown_storage_tag storage_category;
	typedef typename std::size_t size_type;
	typedef typename VType::value_type value_type;
};

/** @internal specialization for unbounded Vectors */
template< typename T >
struct vector_traits< Math::Vector< T, 0 > >
{
	// typedef typename Math::Vector< T, 0 >::base_type base_type;
	typedef typename Math::Vector< T, 0 > self_type;
	typedef dynamic_storage_tag storage_category;
	typedef typename std::size_t size_type;
	typedef T value_type;
	
	// names a function to determine the size of a vector dynamically
	static const size_type size( const self_type& vec ) 
	{
		return vec.size();
	}
};

/** @internal specialization for bounded Vectors */
template< typename T, std::size_t N >
struct vector_traits< Math::Vector< T, N > >
{
	// typedef typename Math::Vector< T, N >::base_type base_type;
	typedef typename Math::Vector< T, N > self_type;
	typedef fixed_storage_tag storage_category;
	typedef typename std::size_t size_type;
	typedef T value_type;

	static const size_type size = N;
};

} } } // namespace Ubitrack::Math::Util

#endif //__UBITRACK_MATH_UTIL_VECTOR_TRAITS_H_INCLUDED__

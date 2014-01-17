

#ifndef __VECTOR_TRAITS_H_INCLUDED__
#define __VECTOR_TRAITS_H_INCLUDED__

//#include <type_traits> //integral_constant <- commented to compile on standard gcc (without c++0x features)
#include <cstddef> // std::size_t
#include "Vector.h"

namespace Ubitrack { namespace Math {

template< class T, T v >
struct constant_value
{
	static const T value = v;
	typedef T value_type;
	typedef constant_value< T, v > type;
	// operator T() { return v; } //<- automatic type conversion, not activated right now
};

typedef constant_value< bool, true > true_type;
typedef constant_value< bool, false > false_type;



template< typename VType >
struct has_fixed_storage
	: public false_type{};
	
template< typename T >
struct has_fixed_storage< Math::Vector< T, 0 > >
	: public false_type{};

template< typename T, std::size_t N >
struct has_fixed_storage< Math::Vector< T, N > >
	: public true_type{};

template< typename VType >
struct has_dynamic_storage
	: public true_type{};

	
template< typename T >
struct has_dynamic_storage< Math::Vector< T, 0 > >
	: public true_type{};

template< typename T, std::size_t N >
struct has_dynamic_storage< Math::Vector< T, N > >
	: public false_type{};


struct fixed_storage_tag {};
struct dynamic_storage_tag {};
struct unknown_storage_tag {};

template< typename VType >
struct vector_traits
{
	typedef VType self_type;
	typedef unknown_storage_tag storage_category;
	typedef typename std::size_t size_type;
	typedef typename VType::value_type value_type;
};



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

}} // namespace Ubitrack::Math

#endif //__VECTOR_TRAITS_H_INCLUDED__

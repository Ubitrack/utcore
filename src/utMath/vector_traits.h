

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


template< typename VType >
struct vector_traits
{
	// typedef typename VType::base_type base_type;
	typedef VType self_type;
	typedef typename std::size_t size_type;
	typedef typename VType::value_type value_type;
};



template< typename T >
struct vector_traits< Math::Vector< T, 0 > >
{
	typedef typename Math::Vector< T, 0 >::base_type base_type;
	typedef typename Math::Vector< T, 0 > self_type;
	typedef typename std::size_t size_type;
	typedef typename T::value_type value_type;
};


template< typename T, std::size_t N >
struct vector_traits< Math::Vector< T, N > >
{
	typedef typename Math::Vector< T, N >::base_type base_type;
	typedef typename Math::Vector< T, N > self_type;
	typedef typename std::size_t size_type;
	typedef T value_type;
	static const size_type size = N;
	// size_type N size_dim;
};


//example usage
// typedef typename Math::vector_traits< Math::Vector< 3, T > >::size_type sizzes;	
// std::cout << "Dynamic " << Math::has_dynamic_storage< Math::Vector< 3, T > >::value << std::endl;


// template< typename T >
// struct vector_traits< CvMat< T > >
// {
	// typedef typename CvMat< T > base_type;
	// typedef typename CvMat< T > self_type;
	// typedef typename int size_type;
	// typedef typename T value_type;
	// // typedef typename elemenstuff< CvMat< T > > element_acess
// };

// struct< typename T >
// elemenstuff
// {
	// T & mat;
	// elementstuff( T & matIn )
	 // : mat( matIn );
	 
	// vec::value_type& operator()( int n )
	// {
		// return mat[ n ];
	// }
	
// }


}} // namespace Ubitrack::Math

#endif //__VECTOR_TRAITS_H_INCLUDED__

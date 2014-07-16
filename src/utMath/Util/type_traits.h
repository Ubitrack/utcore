

#ifndef __UBITRACK_MATH_UTIL_CONSTANT_H_INCLUDED__
#define __UBITRACK_MATH_UTIL_CONSTANT_H_INCLUDED__

//#include <type_traits> //integral_constant <- commented to compile on standard gcc (without c++0x features)

namespace Ubitrack { namespace Math { namespace Util {

/**
 * @internal Constant value struct
 *
 * This functor reassembles the integral constant struct,
 * similar to the c++0x standard and supports the writing
 * of compile time dependent code.
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


/** @internal identifies if a types' size is known at compile time */
template< typename type >
struct has_fixed_storage
	: public false_type{};
	
/** @internal identifies if a types' size is unknown at compile time */
template< typename type >
struct has_dynamic_storage
	: public true_type{};
	
	
struct fixed_storage_tag {};
struct dynamic_storage_tag {};
struct unknown_storage_tag {};


/** 
 * type_trait that defines the precision of a datatypes' underlying built-in type
 *
 * as a default it tries to determine the precision from a defined typedef.
 *
 * it is overwritten for the following built-in-types:
 *   - \c float
 *   - \c double
 *   - \c long \c double
 * where it defaults to the built-in type itself
 */
template< typename T >
struct precision
{
	typedef typename T::value_type type;
};

/// @internal specialization for type float
template<>
struct precision< float >
{
	typedef float type;
};

/// @internal specialization for type double
template<>
struct precision< double >
{
	typedef double type;
};

/// @internal specialization for type long double
template<>
struct precision< long double >
{
	typedef long double type;
};

} } } // namespace Ubitrack::Math::Util

#endif // __UBITRACK_MATH_UTIL_CONSTANT_H_INCLUDED__
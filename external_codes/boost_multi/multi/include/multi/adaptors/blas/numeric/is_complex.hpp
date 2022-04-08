#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020-2021

#ifndef MULTI_ADAPTORS_BLAS_NUMERIC_IS_COMPLEX_HPP
#define MULTI_ADAPTORS_BLAS_NUMERIC_IS_COMPLEX_HPP

#include<complex>
#include<type_traits>

namespace boost::multi::blas::numeric {

using std::true_type;
using std::false_type;

template<class T> auto has_real_fun_aux(T const& t)->decltype(real(t),  true_type{});
                  auto has_real_fun_aux(...       )->decltype(         false_type{});
template<class T> struct has_real_fun : decltype(has_real_fun_aux(std::declval<T>())){};
template<class T> constexpr bool has_real_fun_v = has_real_fun<T>::value;

template<class T> auto has_real_aux(T const& t)->decltype(t.real(),  true_type{});
                  auto has_real_aux(...       )->decltype(          false_type{});
template<class T> struct has_real : decltype(has_real_aux(std::declval<T>())){};
template<class T> constexpr bool has_real_v = has_real<T>::value;

template<class T> auto has_imag_fun_aux(T const& t)->decltype(imag(t),  true_type{});
                  auto has_imag_fun_aux(...       )->decltype(         false_type{});
template<class T> struct has_imag_fun : decltype(has_imag_fun_aux(std::declval<T>())){};
template<class T> constexpr bool has_imag_fun_v = has_imag_fun<T>::value;

template<class T> auto has_imag_aux(T const& t)->decltype(t.imag(),  true_type{});
                  auto has_imag_aux(...       )->decltype(          false_type{});
template<class T> struct has_imag : decltype(has_imag_aux(std::declval<T>())){};
template<class T> constexpr bool has_imag_v = has_imag<T>::value;

template<class T> struct is_complex : std::integral_constant<bool, 
	(has_real_v<T> or has_real_fun_v<T>) and (has_imag_v<T> or has_imag_fun_v<T>)
>{};

template<class V, class T> auto real_is_aux(T const& t)->typename std::is_same<decltype(t.real()), V>;
template<class>            auto real_is_aux(...       )->false_type;
template<class T, class V> struct real_is : decltype(real_is_aux<V>(std::declval<T>())){};

template<class V, class T> auto imag_is_aux(T const& t)->typename std::is_same<decltype(t.imag()), V>;
template<class>            auto imag_is_aux(...       )->false_type;
template<class T, class V> struct imag_is : decltype(imag_is_aux<V>(std::declval<T>())){};

template<class T, class V> struct is_complex_of : std::integral_constant<bool, real_is<T, V>::value and imag_is<T, V>::value>{};

} // end namespace boost::multi::blas::numeric

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS numeric is_complex"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include<thrust/complex.h>

//#include "../../../complex.hpp"
//#include "boost/mpl/list.hpp"

//namespace multi = boost::multi;

//BOOST_AUTO_TEST_CASE(multi_blas_is_complex){
//	namespace numeric = multi::blas::numeric;

//	boost::mpl::for_each<boost::mpl::list<double, float, long double>>([](auto f){
//		using F = decltype(f);
//		static_assert( not numeric::is_complex<F>{}, "!");

//		static_assert( numeric::is_complex<std::complex<F>>{}, "!");
//		static_assert( numeric::is_complex<thrust::complex<F>>{}, "!");
//		static_assert( numeric::is_complex<multi::complex<F>>{}, "!");

//		static_assert( numeric::is_complex_of<std::complex<F>, F>{}, "!");
//		static_assert( not numeric::is_complex_of<F, F>{}, "!");
//	});


//	static_assert( not numeric::is_complex_of<std::complex<double>, float>{}, "!");
//	static_assert( not numeric::is_complex_of<double, float>{}, "!");

//	static_assert( numeric::is_complex<std::complex<double> const&>{}, "!");
//}

//#endif
#endif

// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_NUMERIC_IS_COMPLEX_HPP
#define MULTI_ADAPTORS_BLAS_NUMERIC_IS_COMPLEX_HPP

#include<complex>
#include<type_traits>

namespace boost::multi::blas::numeric {

using std::true_type;
using std::false_type;

template<class T> auto has_real_fun_aux(T const& value) -> decltype(real(value),  true_type{});
                  auto has_real_fun_aux(...           ) -> decltype(             false_type{});
template<class T> struct has_real_fun : decltype(has_real_fun_aux(std::declval<T>())){};
template<class T> constexpr bool has_real_fun_v = has_real_fun<T>::value;

template<class T> auto has_real_aux(T const& value) -> decltype(value.real(),  true_type{});
                  auto has_real_aux(...           ) -> decltype(              false_type{});
template<class T> struct has_real : decltype(has_real_aux(std::declval<T>())){};
template<class T> constexpr bool has_real_v = has_real<T>::value;

template<class T> auto has_imag_fun_aux(T const& value) -> decltype(imag(value),  true_type{});
                  auto has_imag_fun_aux(...           ) -> decltype(         false_type{});
template<class T> struct has_imag_fun : decltype(has_imag_fun_aux(std::declval<T>())){};
template<class T> constexpr bool has_imag_fun_v = has_imag_fun<T>::value;

template<class T> auto has_imag_aux(T const& value) ->decltype(value.imag(),  true_type{});
                  auto has_imag_aux(...           ) ->decltype(          false_type{});
template<class T> struct has_imag : decltype(has_imag_aux(std::declval<T>())){};
template<class T> constexpr bool has_imag_v = has_imag<T>::value;

template<class T> struct is_complex : std::integral_constant<bool, 
	(has_real_v<T> or has_real_fun_v<T>) and (has_imag_v<T> or has_imag_fun_v<T>)
>{};

template<class V, class T> auto real_is_aux(T const& value) -> typename std::is_same<decltype(value.real()), V>;
template<class>            auto real_is_aux(...           ) -> false_type;
template<class T, class V> struct real_is : decltype(real_is_aux<V>(std::declval<T>())){};

template<class V, class T> auto imag_is_aux(T const& value) -> typename std::is_same<decltype(value.imag()), V>;
template<class>            auto imag_is_aux(...           ) -> false_type;
template<class T, class V> struct imag_is : decltype(imag_is_aux<V>(std::declval<T>())){};

template<class T, class V> struct is_complex_of : std::integral_constant<bool, real_is<T, V>::value and imag_is<T, V>::value>{};

} // end namespace boost::multi::blas::numeric

#endif

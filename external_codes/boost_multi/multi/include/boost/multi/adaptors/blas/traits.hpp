// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_TRAITS_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_TRAITS_HPP

#include <complex>
#include <type_traits>  // for enable_if_t, false_type, is_convertible, true...
#include <utility>      // for declval  // IWYU pragma: keep

namespace boost::multi::blas {  // TODO(correaa) include in blas/detail?

// TODO(correaa) : create a BinaryDouble concept?

	template<class F, class=std::enable_if_t<sizeof(F)==sizeof(float ) && std::is_convertible<decltype(std::declval<F&&>()/std::declval<F&&>()), float>{} > >
	auto is_s_aux(F&&) -> std::true_type ;
	auto is_s_aux(...) -> std::false_type;

	template<class T> struct is_s : decltype(is_s_aux(std::declval<T>())) {using archetype = float;};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

	template<class D, class=std::enable_if_t<sizeof(D)==sizeof(double) && std::is_convertible<decltype(std::declval<D&&>()/std::declval<D&&>()), double>{}> >
	auto is_d_aux(D&&) -> std::true_type ;
	auto is_d_aux(...) -> std::false_type;

	template<class T> struct is_d : decltype(is_d_aux(std::declval<T>())) {using archetype = double;};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

	template<class C, class=std::enable_if_t<sizeof(C)==sizeof(std::complex<float>) && is_s<decltype(std::declval<C>().real())>{} && is_s<decltype(std::declval<C>().imag())>{}>>
	auto is_c_aux(C&&) -> std::true_type;
	auto is_c_aux(...) -> std::false_type;

	template<class C> struct is_c : decltype(is_c_aux(std::declval<C>())) {using archetype = std::complex<float>;};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

	template<class Z, class=std::enable_if_t<sizeof(Z)==sizeof(std::complex<double>) && is_d<decltype(std::declval<Z>().real())>{} && is_d<decltype(std::declval<Z>().imag())>{}>>
	auto is_z_aux(Z&&) -> std::true_type ;
	auto is_z_aux(...) -> std::false_type;

	template<class Z> struct is_z : decltype(is_z_aux(std::declval<Z>())) {using archetype = std::complex<double>;};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

} // end namespace boost::multi::blas
#endif

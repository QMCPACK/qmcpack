#ifndef MULTI_ADAPTORS_BLAS_TRAITS_HPP// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
#define MULTI_ADAPTORS_BLAS_TRAITS_HPP
// Copyright 2019-2021 Alfredo A. Correa

#include<complex>
#include<type_traits>

namespace boost::multi::blas{

// TODO(correaa) : create a BinaryDouble concept?

	template<class F, class=std::enable_if_t<sizeof(F)==sizeof(float ) and std::is_convertible<decltype(std::declval<F&&>()/std::declval<F&&>()), float>{}  >>
	auto is_s_aux(F&&) -> std::true_type ;
	auto is_s_aux(...) -> std::false_type;

	template<class T> struct is_s : decltype(is_s_aux(std::declval<T>())){using archetype = float;};

	template<class D, class=std::enable_if_t<sizeof(D)==sizeof(double) and std::is_convertible<decltype(std::declval<D&&>()/std::declval<D&&>()), double>{}>>
	auto is_d_aux(D&&) -> std::true_type ;
	auto is_d_aux(...) -> std::false_type;

	template<class T> struct is_d : decltype(is_d_aux(std::declval<T>())){using archetype = double;};

	template<class C, class=std::enable_if_t<sizeof(C)==sizeof(std::complex<float>) and is_s<decltype(std::declval<C>().real())>{} and is_s<decltype(std::declval<C>().imag())>{}>>
	auto is_c_aux(C&&) -> std::true_type;
	auto is_c_aux(...) -> std::false_type;

	template<class C> struct is_c : decltype(is_c_aux(std::declval<C>())){using archetype = std::complex<float>;};

	template<class Z, class=std::enable_if_t<sizeof(Z)==sizeof(std::complex<double>) and is_d<decltype(std::declval<Z>().real())>{} and is_d<decltype(std::declval<Z>().imag())>{}>>
	auto is_z_aux(Z&&) -> std::true_type ;
	auto is_z_aux(...) -> std::false_type;

	template<class Z> struct is_z : decltype(is_z_aux(std::declval<Z>())){using archetype = std::complex<double>;};

} // end namespace boost::multi::blas
#endif

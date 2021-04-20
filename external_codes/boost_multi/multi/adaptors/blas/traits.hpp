// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_TRAITS_HPP
#define MULTI_ADAPTORS_BLAS_TRAITS_HPP

#include<complex>
#include<type_traits>

namespace boost{
namespace multi::blas{

	template<class F, class=std::enable_if_t<sizeof(F)==sizeof(float ) and std::is_convertible<decltype(std::declval<F&&>()/std::declval<F&&>()), float>{}  >>
	std::true_type  is_s_aux(F&&);
	std::false_type is_s_aux(...);

	template<class T> struct is_s : decltype(is_s_aux(std::declval<T>())){using archetype = float;};

	template<class D, class=std::enable_if_t<sizeof(D)==sizeof(double) and std::is_convertible<decltype(std::declval<D&&>()/std::declval<D&&>()), double>{}>>
	std::true_type  is_d_aux(D&&);
	std::false_type is_d_aux(...);

	template<class T> struct is_d : decltype(is_d_aux(std::declval<T>())){using archetype = double;};

	template<class C, class=std::enable_if_t<sizeof(C)==sizeof(std::complex<float>) and is_s<decltype(std::declval<C>().real())>{} and is_s<decltype(std::declval<C>().imag())>{}>>
	std::true_type  is_c_aux(C&&);
	std::false_type is_c_aux(...);

	template<class C> struct is_c : decltype(is_c_aux(std::declval<C>())){using archetype = std::complex<float>;};

	template<class Z, class=std::enable_if_t<sizeof(Z)==sizeof(std::complex<double>) and is_d<decltype(std::declval<Z>().real())>{} and is_d<decltype(std::declval<Z>().imag())>{}>>
	std::true_type  is_z_aux(Z&&);
	std::false_type is_z_aux(...);

	template<class Z> struct is_z : decltype(is_z_aux(std::declval<Z>())){using archetype = std::complex<double>;};

}

}
#endif


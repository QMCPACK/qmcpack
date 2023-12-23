// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_COMPLEX_TRAITS_HPP
#define MULTI_ADAPTORS_BLAS_COMPLEX_TRAITS_HPP
#pragma once

#include<complex>  // for std::complex
#ifdef __NVCC__
#include<thrust/complex.h>
#endif

namespace boost::multi::blas {

template<class Complex>
struct complex_traits {
	using real_type = typename Complex::real_type;
	constexpr static auto imaginary_unit() { return Complex{real_type{0}, real_type{1}}; }
};

template<class T>
struct complex_traits<std::complex<T>> {
	using real_type = typename ::std::complex<T>::value_type;
	constexpr static auto imaginary_unit() { return ::std::complex<T>{0, 1}; }
};

#ifdef __NVCC__
template<class T>
struct complex_traits<::thrust::complex<T>> {
	using real_type = typename ::thrust::complex<T>::value_type;
	constexpr static auto imaginary_unit() { return std::complex<T>{0, 1}; }
};
#endif


}  // end namespace boost::multi::blas

#endif
// Copyright 2023-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_COMPLEX_TRAITS_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_COMPLEX_TRAITS_HPP

#if defined(__NVCC__) || defined(__HIPCC__)  // defined(__HIP_PLATFORM_AMD__) || defined(__HIP_PLATFORM_NVIDIA__)
#include<thrust/complex.h>
#endif

#include<complex>  // for std::complex

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

#if defined(__NVCC__) || defined(__HIPCC__)  // defined(__HIP_PLATFORM_AMD__) || defined(__HIP_PLATFORM_NVIDIA__)
template<class T>
struct complex_traits<::thrust::complex<T>> {
	using real_type = typename ::thrust::complex<T>::value_type;
	constexpr static auto imaginary_unit() { return std::complex<T>{0, 1}; }
};
#endif

}  // end namespace boost::multi::blas

#endif

// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022-2023 Alfredo A. Correa

#pragma once

#include <multi/array.hpp>

#include<thrust/complex.h>

namespace boost {
namespace multi {

#ifndef NDEBUG
#pragma message "By including this header, the behavior of initialization of thrust::complex<T> in multi::array's changes. ::thrust::complex<T> elements will not be initialized."
#endif

template<class T>
inline constexpr bool force_element_trivial_default_construction<::thrust::complex<T>> = std::is_trivially_default_constructible_v<T>;

template<class T> struct is_trivially_default_constructible<::thrust::complex<T>> : std::is_trivially_default_constructible<T> {};

static_assert(is_trivially_default_constructible<::thrust::complex<double>>::value);

}
}

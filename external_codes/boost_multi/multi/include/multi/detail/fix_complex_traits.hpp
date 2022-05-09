// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

#pragma once

#include "../detail/type_traits.hpp"

#include<complex>

namespace boost {
namespace multi {

#pragma message "By including this header, the behavior of initialization of std::complex<T> in multi::array's changes. std::complex<T>  elements will not be initialized."

template<class T> struct is_trivially_default_constructible<std::complex<T>> : std::is_trivially_default_constructible<T> {};

}
}

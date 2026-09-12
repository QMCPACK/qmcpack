// Copyright 2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_FILLING_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_FILLING_HPP
#pragma once

#include "boost/multi/detail/config/UNREACHABLE.hpp"  // for BOOST_MULTI_UNREACHABLE

// TODO(correaa)  #include "multi/blas/filling.hpp"

namespace boost::multi::lapack {

enum class filling : char {
	lower = 'U',
	upper = 'L',
};

#ifdef __NVCC__  // in place of global -Xcudafe "--diag_suppress=implicit_return_from_non_void_function"
#pragma nv_diagnostic push
#pragma nv_diag_suppress = implicit_return_from_non_void_function  // nvcc EDG front end doesn't see BOOST_MULTI_UNREACHABLE as noreturn in MSVC-host mode
#endif
#ifdef __NVCOMPILER
#pragma diagnostic push
#pragma diag_suppress = implicit_return_from_non_void_function
#endif

inline auto flip(filling side) -> filling {
	switch(side) {
	case filling::lower: return filling::upper;
	case filling::upper: return filling::lower;
	}
	BOOST_MULTI_UNREACHABLE();  // LCOV_EXCL_LINE
}

#ifdef __NVCOMPILER
#pragma diagnostic pop
#endif
#ifdef __NVCC__
#pragma nv_diagnostic pop
#endif

inline auto operator-(filling side) -> filling { return flip(side); }
inline auto operator+(filling side) -> filling { return side; }

}  // namespace boost::multi::lapack

#endif

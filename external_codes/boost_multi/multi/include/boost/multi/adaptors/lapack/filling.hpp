// Copyright 2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_FILLING_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_FILLING_HPP
#pragma once

// TODO(correaa)  #include "multi/blas/filling.hpp"

namespace boost::multi::lapack {

enum class filling : char {
	lower = 'U',
	upper = 'L',
};

inline auto flip(filling side) -> filling {
	switch(side) {
	case filling::lower: return filling::upper;
	case filling::upper: return filling::lower;
	}
	__builtin_unreachable();  // LCOV_EXCL_LINE
}

inline auto operator-(filling side) -> filling { return flip(side); }
inline auto operator+(filling side) -> filling { return side; }

}  // namespace boost::multi::lapack

#endif

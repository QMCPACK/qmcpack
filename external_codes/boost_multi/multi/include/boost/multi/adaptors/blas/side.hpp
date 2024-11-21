// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_SIDE_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_SIDE_HPP
#pragma once

namespace boost::multi::blas {

enum class side : char {
	left  = 'L',
	right = 'R'
};

inline auto swap(side sid) noexcept -> side {
	switch(sid) {
		case side::left : return side::right;
		case side::right: return side::left ;
	}  // __builtin_unreachable();  // LCOV_EXCL_LINE
	return {};
}

} // end namespace boost::multi::blas
#endif

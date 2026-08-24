// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_SIDE_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_SIDE_HPP

namespace boost::multi::blas {

enum class side : char {
	left  = 'L',
	right = 'R'
};

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wswitch-default"
#endif

inline auto swap(side sid) noexcept -> side {
	switch(sid) {  // NOLINT(clang-diagnostic-switch-default)
		case side::left : return side::right;
		case side::right: return side::left ;
	}  // __builtin_unreachable();  // LCOV_EXCL_LINE
	return side::left;
	// return {};
}

inline auto iden(side sid) noexcept -> side {
	return sid;
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif


} // end namespace boost::multi::blas
#endif

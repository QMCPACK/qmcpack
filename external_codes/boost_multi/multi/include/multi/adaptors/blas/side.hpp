// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_SIDE_HPP
#define MULTI_ADAPTORS_BLAS_SIDE_HPP

namespace boost::multi::blas {

enum side : char {
	left  = 'L',
	right = 'R'//,
//  pre_multiply = 'R',
//  post_multiply = 'L'
};

inline auto swap(side sid) -> side {
	switch(sid) {
		case side::left : return side::right;
		case side::right: return side::left ;
	} __builtin_unreachable();  // LCOV_EXCL_LINE
}

} // end namespace boost::multi::blas

#endif

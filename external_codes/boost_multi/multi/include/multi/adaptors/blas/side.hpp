#ifndef MULTI_ADAPTORS_BLAS_SIDE_HPP // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
#define MULTI_ADAPTORS_BLAS_SIDE_HPP
// Copyright 2019-2021 Alfredo A. Correa

namespace boost::multi::blas {

enum side : char {
	left  = 'L', 
	right = 'R'//,
//	pre_multiply = 'R', 
//	post_multiply = 'L'
};

inline auto swap(side sid) -> side {
	switch(sid) {
		case side::left : return side::right;
		case side::right: return side::left ;
	} __builtin_unreachable();
}

} // end namespace boost::multi::blas

#endif

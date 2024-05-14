// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_SWAP_HPP
#define MULTI_ADAPTORS_BLAS_SWAP_HPP
#pragma once

#include <multi/adaptors/blas/core.hpp>

namespace boost::multi::blas {

using core::swap;

template<class It1, class Size, class It2>
auto swap_n(It1 first, Size count, It2 first2) -> It2 {
	blas::default_context_of(base(first))->swap(count, base(first), stride(first), base(first2), stride(first2));
	return first2 + count;
}

template<class It1, class It2>
auto swap(It1 first, It2 last, It2 first2) -> It2 {
	assert(stride(first) == stride(last));
	return swap_n(first, last - first, first2);
}

template<class X1D, class Y1D>
auto swap(X1D&& x, Y1D&& y) -> Y1D&& {  // NOLINT(readability-identifier-length) x, y conventional blas names
	assert( size(x) == size(y) );
//  assert( offset(x) == 0 and offset(y) == 0 );
	swap( begin(x), end(x), begin(y) );
	return std::forward<Y1D>(y);
}

template<class X1D, class Y1D>
auto swap(X1D const& x, Y1D const& y) = delete;  // NOLINT(readability-identifier-length) x, y conventional blas names

template<class X1D, class Y1D>
auto operator^(X1D&& x, Y1D&& y) {  // NOLINT(readability-identifier-length) BLAS naming
	blas::swap(x, y);
	return std::tie(x, y);
}

namespace operators {
	using blas::operator^;
}  // end namespace operators

} // end namespace boost::multi::blas
#endif

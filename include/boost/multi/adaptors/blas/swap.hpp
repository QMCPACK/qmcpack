// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_SWAP_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_SWAP_HPP
#pragma once

#include "boost/multi/adaptors/blas/core.hpp"

namespace boost::multi::blas {

using core::swap;

template<class It1, class Size, class It2>
auto swap_n(It1 first, Size count, It2 first2) -> It2 {
	blas::default_context_of(base(first))->swap(count, base(first), stride(first), base(first2), stride(first2));
	return first2 + count;
}

template<class It1, class It2>
auto swap(It1 first, It2 last, It2 first2) noexcept -> It2 {
	assert(stride(first) == stride(last));
	return swap_n(first, last - first, first2);
}

template<class X1D, class Y1D>
auto swap(X1D&& x, Y1D&& y) noexcept(false) -> Y1D&& {  // NOLINT(readability-identifier-length) x, y conventional blas names, // NOSONAR(cpp:S5018) this swap can "fail" if sizes do not match
	assert( size(x) == size(y) );
	swap( std::begin(x), std::end(std::forward<X1D>(x)), std::begin(y) );
	return std::forward<Y1D>(y);
}

template<class X1D, class Y1D>
auto swap(X1D const&, Y1D const&) noexcept(false) = delete;  // NOSONAR(cpp:S5018) this swap can "fail" if sizes do not match

template<class X1D, class Y1D>
auto operator^(X1D&& x, Y1D&& y) {  // NOLINT(readability-identifier-length) BLAS naming
	blas::swap(x, y);
	return std::tie(std::forward<X1D>(x), std::forward<Y1D>(y));  // or use std::forward_as_tuple ?
}

namespace operators {
	using blas::operator^;
}  // end namespace operators

} // end namespace boost::multi::blas
#endif

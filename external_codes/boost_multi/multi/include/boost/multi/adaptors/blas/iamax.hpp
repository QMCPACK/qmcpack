// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_IAMAX_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_IAMAX_HPP

#include <boost/multi/adaptors/blas/core.hpp>

namespace boost::multi::blas {

template<class It, class Size>
auto iamax_n(It first, Size n) {
	using core::iamax;
	return iamax(n, base(first), stride(first));  // if you get an error here make sure that you are including (and linking) the appropriate BLAS backend for your memory type
}

template<class It>
auto iamax(It first, It last)
	-> decltype(iamax_n(first, std::distance(first, last))) {
	return iamax_n(first, std::distance(first, last));
}

template<class X1D>
auto iamax(X1D const& x)  // NOLINT(readability-identifier-length) x conventional blas name
	-> decltype(iamax(begin(x), end(x))) {
	assert(! offset(x));
	return iamax(begin(x), end(x));
}

template<class X1D>
auto amax(X1D const& x) {  // NOLINT(readability-identifier-length) x conventional blas name
	return begin(x) + iamax(x);
}

}  // end namespace boost::multi::blas

#endif

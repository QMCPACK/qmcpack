// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_SCAL_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_SCAL_HPP
#pragma once

#include <boost/multi/adaptors/blas/core.hpp>
// IWYU pragma: no_include "boost/multi/adaptors/blas/traits.hpp"  // for blas, multi

// IWYU pragma: no_include <type_traits>  // for declval  // needed by iwyu-clang-linux
#include <utility>  // for forward, declval

namespace boost::multi::blas {

using core::scal;

template<class It, class Size>
auto scal_n(typename It::element a, It first, Size count) {  // NOLINT(readability-identifier-length) conventional BLAS naming
	auto ctxt = blas::default_context_of(first.base());
	ctxt->scal(count, &a, first.base(), first.stride());
}

template<class Scalar, class It1D>
auto scal(Scalar const& a, It1D first, It1D last)  // NOLINT(readability-identifier-length) conventional BLAS naming
->decltype(blas::scal_n(a, first, last - first)) {  // NOLINT(fuchsia-default-arguments-calls) allow a possible double -> complex conversion (with default 0 imag part)
	return blas::scal_n(a, first, last - first); }  // NOLINT(fuchsia-default-arguments-calls) same

template<class Scalar, class X1D>  // don't do this: ", typename Elem = typename X1D::element_type>"
auto scal(Scalar const& a, X1D&& x)  // NOLINT(readability-identifier-length) conventional BLAS naming
->decltype(blas::scal(a, x.begin(), x.end()), std::forward<X1D>(x)) {
	return blas::scal(a, x.begin(), x.end()), std::forward<X1D>(x); }

template<class A>
class scal_range {
	A alpha_;

 public:
	using scalar_type = A;
	explicit scal_range(A const& alpha) : alpha_{alpha} {}
	template<class X1D>
	friend auto operator*=(X1D&& x, scal_range const& self)  // NOLINT(readability-identifier-length) conventional BLAS naming
	->decltype(std::forward<X1D>(scal(std::declval<scalar_type const&>(), std::forward<X1D>(x)))) {
		return std::forward<X1D>(scal(self.alpha_, std::forward<X1D>(x)));}
};

template<class A> auto scal(A const& array) {return scal_range<A>{array};}

namespace operators {
	template<class Scalar, class X>
	auto operator*=(X&& x, Scalar const& alpha) -> X&& {  // NOLINT(readability-identifier-length) conventional BLAS naming
		return blas::scal(alpha, std::forward<X>(x));
	}
}  // end namespace operators

}  // end namespace boost::multi::blas
#endif

#ifndef MULTI_ADAPTORS_BLAS_SCAL_HPP  // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
#define MULTI_ADAPTORS_BLAS_SCAL_HPP
// Copyright 2019-2021 Alfredo A. Correa

#include "../blas/core.hpp"

namespace boost::multi::blas{

using core::scal;

template<class Scalar, class It, class Size>
auto scal_n(Scalar const& a, It first, Size count)  // NOLINT(readability-identifier-length) conventional BLAS naming
->decltype(scal(count, &a, first.base(), first.stride()), void()) {
	       scal(count, &a, first.base(), first.stride());         }

template<class Scalar, class It1D>
auto scal(Scalar const& a, It1D first, It1D last)  // NOLINT(readability-identifier-length) conventional BLAS naming
->decltype(blas::scal_n(a, first, last - first)){
	return blas::scal_n(a, first, last - first);}

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
	->decltype(std::forward<X1D>(scal(std::declval<scalar_type const&>(), x))) {
		return std::forward<X1D>(scal(self.alpha_, x));}
};

template<class A> auto scal(A const& array) {return scal_range<A>{array};}

} // end namespace boost::multi::blas

#endif

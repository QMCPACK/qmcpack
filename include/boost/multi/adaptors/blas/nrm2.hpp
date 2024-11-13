// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_NRM2_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_NRM2_HPP
#pragma once

#include <boost/multi/adaptors/blas/core.hpp>

#include <boost/multi/array.hpp>

#include<complex>  // std::norm

namespace boost::multi::blas {

using core::nrm2;

using multi::base;
using std::norm;  // nvcc11 needs using std::FUNCTION and the FUNCTION (and it works in clang, gcc, culang, icc)

template<class Context, class XIt, class Size, class RPtr>
auto nrm2_n(Context&& ctxt, XIt x_first, Size count, RPtr rp) {
	std::forward<Context>(ctxt)->nrm2(count, x_first.base(), x_first.stride(), rp);
}

template<class It, class Size, class A0D>
auto nrm2_n(It const& x, Size n, A0D res)  // NOLINT(readability-identifier-length) conventional BLAS naming
//->decltype(blas::default_context_of(x.base())->nrm2(n, x.base(), x.stride(), res), std::next(res)) {  // NOLINT(fuchsia-default-arguments-calls)
{   return blas::default_context_of(x.base())->nrm2(n, x.base(), x.stride(), res), std::next(res); }  // NOLINT(fuchsia-default-arguments-calls)

template<class Context, class X1D, class R,
	std::enable_if_t<! multi::has_base<std::decay_t<R>>::value, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto nrm2(Context ctxt, X1D const& x, R&& res) -> R&& {  // NOLINT(readability-identifier-length) res = \sum_i x_i y_i
	return blas::nrm2_n(ctxt, x.begin(), size(x), &res), std::forward<R>(res);
}

template<class Context, class X1D, class R,
	std::enable_if_t<multi::has_base<std::decay_t<R>>::value, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto nrm2(Context ctxt, X1D const& x, R&& res) -> R&& {  // NOLINT(readability-identifier-length) res = \sum_i x_i y_i
	return blas::nrm2_n(ctxt, begin(x), size(x), res.base()), std::forward<R>(res);
}

template<class X1D, class R>
auto nrm2(X1D const& x, R&& res) -> R&& {  // NOLINT(readability-identifier-length) BLAS naming
	auto ctxtp = blas::default_context_of(x.base());
	return blas::nrm2(ctxtp, x, std::forward<R>(res));
}

template<class ItX, class Size>
class nrm2_ptr {
	ItX  x_first_;
	Size count_;

 protected:
	nrm2_ptr(ItX x_first, Size count) : x_first_{x_first}, count_{count} {}

 public:
	explicit operator bool() const {return true;}

	template<class ItOut, class Size2>
	friend constexpr auto copy_n(nrm2_ptr first, Size2 count, ItOut d_first) {
//  ->decltype(blas::nrm2_n(std::declval<ItX>(), Size2{}     , d_first), d_first + count) {
		assert(count == 1); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		return blas::nrm2_n(first.x_first_     , first.count_, d_first), d_first + count; }

	template<class ItOut, class Size2>
	friend constexpr auto uninitialized_copy_n(nrm2_ptr first, Size2 count, ItOut d_first)
	->decltype(blas::nrm2_n(std::declval<ItX>(), Size2{}     , d_first), d_first + count) {assert(count == 1);
		return blas::nrm2_n(first.x_first_     , first.count_, d_first), d_first + count; }

	template<class ItOut, class Size2>
	static constexpr auto uninitialized_copy_n(nrm2_ptr first, Size2 count, ItOut d_first)
	->decltype(blas::nrm2_n(std::declval<ItX>(), Size2{}     , d_first), d_first + count) {assert(count == 1);
		return blas::nrm2_n(first.x_first_     , first.count_, d_first), d_first + count; }
};

template<class X, class Ptr = nrm2_ptr<typename X::const_iterator, typename X::size_type>>
struct nrm2_ref : private Ptr {
	using decay_type = decltype(norm(std::declval<typename X::value_type>()));
	explicit nrm2_ref(X const& x) : Ptr{begin(x), size(x)} {}  // NOLINT(readability-identifier-length) BLAS naming

	constexpr auto operator&() const& -> Ptr const& {return *this;}  // NOLINT(google-runtime-operator) reference type  //NOSONAR

	auto decay() const -> decay_type {decay_type ret; copy_n(operator&(), 1, &ret); return ret;}  // NOLINT(fuchsia-default-arguments-calls) complex
	operator decay_type()       const { return decay(); }  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,hicpp-explicit-conversion) //NOSONAR to allow terse syntax
// #if ! defined(__CUDACC__) || ! defined(__INTEL_COMPILER)
//  friend auto operator*(decay_type const& lhs, dot_ref const& self) {return lhs*self.decay();}
// #endif
	auto operator+() const -> decay_type { return decay(); }

	// auto operator==(nrm2_ref const& other) const { return decay() == other.decay(); }
	// auto operator!=(nrm2_ref const& other) const { return decay() != other.decay(); }

	friend auto operator==(nrm2_ref const& self, nrm2_ref const& other) { return self.decay() == other.decay(); }
	friend auto operator!=(nrm2_ref const& self, nrm2_ref const& other) { return self.decay() != other.decay(); }
};

template<class X>
[[nodiscard]]
auto nrm2(X const& x) {  // NOLINT(readability-identifier-length) BLAS naming
	return nrm2_ref<X>{x};
}

namespace operators {
	using std::norm;
	template<class A1D, class Real = decltype(norm(std::declval<typename A1D::value_type>()))>//decltype(norm(std::declval<typename A1D::value_type>()))>
	[[nodiscard]] auto operator^(A1D const& array, int n)
	->decltype(std::pow(Real{blas::nrm2(array)}, n)) {
		return std::pow(Real{blas::nrm2(array)}, n); }

	template<class A1D>
	[[nodiscard]] auto abs(A1D const& array) {
		return blas::nrm2(array);
	}

	template<class A1D>
	[[nodiscard]] auto norm(A1D const& array) {
		auto const sqrt = +blas::nrm2(array);
		return sqrt*sqrt;
	}

} // end namespace operators

} // end namespace boost::multi::blas

#endif

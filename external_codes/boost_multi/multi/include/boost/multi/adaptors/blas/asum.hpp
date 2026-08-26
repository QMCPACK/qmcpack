// Copyright 2019-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_ASUM_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_ASUM_HPP
#pragma once

#include "boost/multi/adaptors/blas/core.hpp"

namespace boost::multi::blas {

template<class It, typename Size, class A0D>
auto asum_n(It first, Size n, A0D res)
->decltype(blas::default_context_of(base(first))->asum(n, base(first), stride(first), res), std::next(res)) {
	return blas::default_context_of(base(first))->asum(n, base(first), stride(first), res), std::next(res); }

using std::begin; using std::end;

template<class X1D, class A0D>
auto asum(X1D const& x, A0D&& res)  // NOLINT(readability-identifier-length) x conventional blas name
//->decltype(asum_n(x.begin(), x.size(), &res)) {
{   return asum_n(std::begin(x), x.size(), &std::forward<A0D>(res)); }

template<class A1D>
struct asum_ptr {
	A1D const* xp_;  // NOLINT(misc-non-private-member-variables-in-classes)

	explicit operator bool() const {return xp_;}

	template<class ItOut, class Size2>
	friend auto copy_n(asum_ptr first, [[maybe_unused]] Size2 count, ItOut d_first)
	->decltype(blas::asum_n(typename A1D::iterator{}, typename A1D::size_type{}, d_first)) {assert(count == 1);
		return blas::asum_n(first.xp_->begin()      , first.xp_->size()        , d_first); }

	template<class... As>
	friend auto uninitialized_copy_n(asum_ptr first, As... as) {return copy_n(first, as...);}

	template<class... As>
	static auto uninitialized_copy_n(asum_ptr first, As... as) {return copy_n(first, as...);}
};

template<class A1D>
[[nodiscard]]
auto asum(A1D const& x) {  // NOLINT(readability-identifier-length) BLAS naming
	struct ref {
		A1D const& x_;  // NOLINT(misc-non-private-member-variables-in-classes,cppcoreguidelines-avoid-const-or-ref-data-members)
		auto operator&() const& {return asum_ptr<A1D>{&x_};}  // NOLINT(google-runtime-operator) reference type //NOSONAR
		using decay_type = decltype(abs(std::declval<typename A1D::value_type>()));
		operator decay_type() const {decay_type ret; blas::asum(x_, ret); return ret;}  //NOSONAR // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)  allow terse syntax double a = asum(v);
		auto operator+() const -> decay_type {return operator decay_type();}
	};

	return ref{x};
}

namespace operators {
	static constexpr double threshold = 1.0e-12;

	using zero_type = void*****;

	template<class A1D>
	auto operator==(A1D const& self, [[maybe_unused]] zero_type zero) -> bool {
		assert( zero == nullptr );
		return blas::asum(self) < threshold;
	}

	template<class A1D>
	auto operator!=(A1D const& self, [[maybe_unused]] zero_type zero) -> bool {
		assert( zero == nullptr );
		return blas::asum(self) > threshold;
	}

	template<class A1D>
	auto contains_nan(A1D const& self) -> bool {
		return std::isnan(blas::asum(self));
	}

	template<class A1D>
	auto isnan(A1D const& self) -> bool {
		return contains_nan(self);
	}

	template<class A1D>
	auto isfinite(A1D const& self) -> bool {
		return std::isfinite(blas::asum(self));
	}

	template<class A1D>
	auto isinf(A1D const& self) -> bool {
		return std::isinf(blas::asum(self));
	}
}  // end namespace operators

}  // end namespace boost::multi::blas

#endif

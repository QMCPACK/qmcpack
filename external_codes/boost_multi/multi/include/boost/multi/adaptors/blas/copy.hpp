// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_COPY_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_COPY_HPP
#pragma once

#include <boost/multi/adaptors/blas/core.hpp>  // for copy, default_context_of  // IWYU pragma: export
// IWYU pragma: no_include "boost/multi/adaptors/blas/core.hpp"  // bug in iwyu 18.1.8?

#include <cassert>   // for assert
#include <iterator>  // for iterator_traits, outpu...
#include <utility>   // for forward

namespace boost::multi::blas {

using core::copy;

template<class It, typename Size, class OutIt>
auto copy_n(It first, Size n, OutIt d_first)
	-> decltype(blas::default_context_of(first.base())->copy(n, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n) {
	/**/ return blas::default_context_of(first.base())->copy(n, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n;
}

template<class X1D, class Y1D>
auto copy(X1D const& x, Y1D&& y)  // NOLINT(readability-identifier-length) BLAS naming
	-> decltype(blas::copy_n(x.begin(), size(x), y.begin()), std::forward<Y1D>(y)) {
	assert((x.size() == y.size()));  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : assert
	/**/ return blas::copy_n(x.begin(), x.size(), y.begin()), std::forward<Y1D>(y);
}

template<class It>
struct copy_it {
	It it_;  // NOLINT(misc-non-private-member-variables-in-classes)

	using difference_type   = typename std::iterator_traits<It>::difference_type;
	using value_type        = typename std::iterator_traits<It>::value_type;
	using pointer           = void;
	using reference         = void;
	using iterator_category = std::output_iterator_tag;
	using iterator_type     = copy_it;

	friend auto operator-(copy_it const& c1, copy_it const& c2) { return c1.it_ - c2.it_; }

	template<class It1DOut>
	friend constexpr auto copy_n(copy_it first, difference_type count, It1DOut result) -> It1DOut {
		return blas::copy_n(first.it_, count, result);
	}
	template<class It1DOut>
	friend constexpr auto copy(copy_it first, copy_it last, It1DOut d_first) -> It1DOut {
		return copy_n(first, distance(first, last), d_first);
	}
	template<class It1DOut>
	friend constexpr auto uninitialized_copy(copy_it first, copy_it last, It1DOut d_first) -> It1DOut {
		return copy_n(first, distance(first, last), d_first);
	}
	friend constexpr auto distance(copy_it const& self, copy_it const& other) -> difference_type {
		return other.it_ - self.it_;
	}
	constexpr auto operator*() const -> value_type { return *it_; }
};

template<class A1D> [[nodiscard]]
auto copy(A1D const& x) {  // NOLINT(readability-identifier-length) BLAS naming
	struct ref {
		A1D const& x_;  // NOLINT(misc-non-private-member-variables-in-classes,cppcoreguidelines-avoid-const-or-ref-data-members)
		using iterator = copy_it<typename A1D::const_iterator>;
		auto begin() const { return iterator{x_.begin()}; }
		auto end() const { return iterator{x_.end()}; }
		auto size() const { return x_.size(); }
		auto extensions() const { return x_.extensions(); }
	};
	return ref{x};
}

namespace operators {

template<class A1D, class B1D>
auto operator<<(A1D&& lhs, B1D const& rhs) -> A1D&& {
	return boost::multi::blas::copy(rhs, std::forward<A1D>(lhs));
}

}  // end namespace operators

}  // end namespace boost::multi::blas

#endif

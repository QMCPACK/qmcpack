// Copyright 2020-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_COPY_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_COPY_HPP

#include "boost/multi/adaptors/blas/core.hpp"  // IWYU pragma: keep

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
	assert((x.size() == y.size()));
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

	explicit copy_it(It it) : it_{it} {}

	copy_it() = default;

	copy_it(copy_it const&)     = default;
	copy_it(copy_it&&) noexcept = default;

	auto operator=(copy_it const&) -> copy_it&     = default;
	auto operator=(copy_it&&) noexcept -> copy_it& = default;

	~copy_it() = default;

	friend auto operator-(copy_it const& c1, copy_it const& c2) { return c1.it_ - c2.it_; }

	auto operator==(copy_it const& other) const -> bool { return it_ == other.it_; }
	auto operator!=(copy_it const& other) const -> bool { return it_ != other.it_; }

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

	constexpr auto operator++() -> copy_it&;  // void {};  // -> copy_it& { assert(0); /*++it_;*/ return *this; }
	constexpr auto operator--() -> copy_it&;  // void {};  // -> copy_it& { assert(0); /*--it_;*/ return *this; }

	constexpr auto operator++(int) -> copy_it;
	constexpr auto operator--(int) -> copy_it;
};

template<class A1D> [[nodiscard]]
auto copy(A1D const& x) {  // NOLINT(readability-identifier-length) BLAS naming
	struct ref {
		A1D const& x_;  // NOLINT(misc-non-private-member-variables-in-classes,cppcoreguidelines-avoid-const-or-ref-data-members)
		using iterator = copy_it<typename A1D::const_iterator>;
		auto begin() const { return iterator{x_.begin()}; }
		auto end() const { return iterator{x_.end()}; }
		auto size() const { return x_.size(); }
		auto extensions() const { return x_.extents(); }
		[[nodiscard]] constexpr auto extents() const { return x_.extents(); }
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

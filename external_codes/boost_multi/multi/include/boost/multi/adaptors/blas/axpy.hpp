// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_AXPY_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_AXPY_HPP
#pragma once

#include <boost/multi/adaptors/blas/core.hpp>
// #include <boost/multi/adaptors/complex.hpp>

#include <boost/multi/array_ref.hpp>

#include <cassert>                            // for assert
#include <complex>                             // for complex
#include <iterator>                            // for iterator_traits, rando...
#include <type_traits>                         // for enable_if_t, decay_t
#include <utility>                             // for forward, declval

#define BOOST_MULTI_DECLRETURN(ExpR) -> decltype(ExpR) {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing
#define BOOST_MULTI_JUSTRETURN(ExpR)                   {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

namespace boost::multi::blas {

using core::axpy;

template<class It1, class Size, class OutIt>
auto axpy_n(typename It1::value_type alpha, It1 first, Size n, OutIt d_first)
->decltype(axpy(n, &alpha, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n) {
	return axpy(n, &alpha, first.base(), first.stride(), d_first.base() , stride(d_first) ), d_first + n; }

template<class Context, class It1, class Size, class OutIt>//, class=std::enable_if_t<is_context<decltype(*Context{})>{}>>
auto axpy_n(Context ctxt, typename It1::value_type alpha, It1 first, Size n, OutIt d_first)
->decltype(ctxt->axpy(n, &alpha, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n) {
	return ctxt->axpy(n, &alpha, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n; }

template<class Context, class X1DIt, class Y1D, typename = decltype( std::declval<Y1D&&>()[0] = 0.0, *X1DIt{} )>
auto axpy(Context ctxt, typename X1DIt::element alpha, X1DIt x, Y1D&& y)  // NOLINT(readability-identifier-length) conventional BLAS names
->decltype(                                                 std::forward<Y1D>(y)) {
	return axpy_n(ctxt, alpha,   x,   size(y),   begin(y)), std::forward<Y1D>(y); }

template<class Context, class X1D, class Y1D, typename = decltype( std::declval<Y1D&&>()[0] = 0.0, size(std::declval<X1D const&>()) )>
auto axpy(Context ctxt, typename X1D::element alpha, X1D const& x, Y1D&& y)  // NOLINT(readability-identifier-length) conventional BLAS names
->decltype(                                                        std::forward<Y1D>(y)) { assert(x.size() == y.size() );
	return axpy_n(ctxt, alpha,   x.begin(),   y.size(),   y.begin()), std::forward<Y1D>(y); }

template<class X1D, class Y1D, typename = decltype( std::declval<Y1D&&>()[0] = 0.0 )>
auto axpy(typename X1D::element alpha, X1D const& x, Y1D&& y)  // NOLINT(readability-identifier-length) conventional BLAS names
-> decltype(auto)
{
	auto ctxtp = blas::default_context_of(x.base());
	return boost::multi::blas::axpy(ctxtp, alpha, x, std::forward<Y1D>(y));
}

template<class X1D, class Y1D>
auto axpy(X1D const& x, Y1D&& y) -> Y1D&& {  // NOLINT(readability-identifier-length) conventional BLAS names
	return axpy(+1.0, x, std::forward<Y1D>(y));
}

template<class Context, class X1D, class Y1D, std::enable_if_t<is_context<Context>{}> >
auto axpy(Context&& ctxt, X1D const& x, Y1D&& y) -> Y1D&& {  // NOLINT(readability-identifier-length) conventional BLAS names
	return axpy(std::forward<Context>(ctxt), +1.0, x, std::forward<Y1D>(y));
}

template<class Context, class Scale, class ItX>
class axpy_iterator {
	Context ctxt_;
	Scale alpha_;
	ItX x_begin_;

 public:
	axpy_iterator(Context ctxt, Scale alpha, ItX x_begin)
	: ctxt_{ctxt}, alpha_{alpha}, x_begin_{x_begin} {}

	using difference_type = typename std::iterator_traits<ItX>::difference_type;
	using value_type = typename std::iterator_traits<ItX>::value_type;
	using pointer = void;
	using reference = void;
	using iterator_category = std::random_access_iterator_tag;

	friend auto operator-(axpy_iterator const& self, axpy_iterator const& other) -> difference_type {
		assert(self.alpha_ == other.alpha_);
		return self.x_begin_ - other.x_begin_;
	}

	template<class It1DOut>
	friend auto copy_n(axpy_iterator first, difference_type count, It1DOut result) {
		blas::axpy_n(first.ctxt_, first.alpha_, first.x_begin_, count, result);  // NOLINT(fuchsia-default-arguments-calls)
		return result + count;
	}
	template<class It1DOut>
	friend auto copy(axpy_iterator first, axpy_iterator last, It1DOut result){return copy_n(first, last - first, result);}
};

template<class Context, class Scale, class ItX>
class axpy_range {
	Context ctxt_;
	Scale alpha_;
	ItX x_begin_;
	size_type count_;

 public:
	axpy_range(axpy_range const&) = delete;
	axpy_range(axpy_range&&) noexcept = delete;
	~axpy_range() = default;
	auto operator=(axpy_range const&) -> axpy_range& = delete;
	auto operator=(axpy_range&&) noexcept -> axpy_range& = delete;

	axpy_range(Context ctxt, Scale alpha, ItX x_first, ItX x_last)
	: ctxt_{ctxt}, alpha_{alpha}, x_begin_{x_first}, count_{x_last - x_first} {}

	using iterator = axpy_iterator<Context, Scale, ItX>;

	auto begin() const -> iterator{ return {ctxt_, alpha_, x_begin_         }; }
	auto end()   const -> iterator{ return {ctxt_, alpha_, x_begin_ + count_}; }

	auto size() const -> size_type { return end() - begin(); }
	auto extensions() const { return extensions_t<1>{ {0, size()} }; }

	template<class Other>
	friend auto operator+=(Other&& other, axpy_range const& self) -> Other&& {
		assert(other.size() == self.count_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : bug in clang-tidy https://reviews.llvm.org/D31130
		blas::axpy_n(self.ctxt_, +static_cast<typename ItX::value_type>(self.alpha_), self.x_begin_, self.count_, other.begin());  // NOLINT(fuchsia-default-arguments-calls)
		return std::forward<Other>(other);
	}
	template<class Other>
	friend auto operator-=(Other&& other, axpy_range const& self) -> Other&& {
		assert(other.size() == self.count_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : bug in clang-tidy https://reviews.llvm.org/D31130
		blas::axpy_n(self.ctxt_, -static_cast<typename ItX::value_type>(self.alpha_), self.x_begin_, self.count_, other.begin());  // NOLINT(fuchsia-default-arguments-calls)
		return std::forward<Other>(other);
	}
	auto operator*=(Scale s) & -> axpy_range& {alpha_ *= s; return *this;}  // NOLINT(readability-identifier-length) conventional BLAS naming
};

template<class Context, class Scalar, class X1D,
	class=std::enable_if_t<is_context<Context>{}>>  // NOLINT(modernize-use-constraints) for C++20
auto axpy(Context&& ctxt, Scalar a, X1D const& x)  // NOLINT(readability-identifier-length) conventional BLAS naming
-> axpy_range<Context, Scalar, typename X1D::const_iterator> {  // NOLINT(readability-identifier-length) conventional BLAS naming
	return {std::forward<Context>(ctxt), a, begin(x), end(x)};
}

template<class Scalar, class X1D>
auto axpy(Scalar a, X1D const& x)  // NOLINT(readability-identifier-length) conventional BLAS naming
{
	auto ctxtp = blas::default_context_of(x.base());
	return axpy_range<decltype(ctxtp), Scalar, typename X1D::const_iterator>{ctxtp, a, begin(x), end(x)};
}

template<class AA, class X>
class scaled {
	AA a_;
	X const& x_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

 public:
	scaled(AA a, X const& x) : a_{a}, x_{x} {}  // NOLINT(readability-identifier-length) conventional BLAS naming
	template<class Y1D>
	friend auto operator+=(Y1D&& y, scaled const& ax) {return axpy(+ax.a_, ax.x_, std::forward<Y1D>(y));}  // NOLINT(readability-identifier-length) conventional BLAS naming
	template<class Y1D>
	friend auto operator-=(Y1D&& y, scaled const& ax) {return axpy(-ax.a_, ax.x_, std::forward<Y1D>(y));}  // NOLINT(readability-identifier-length) conventional BLAS naming
};

namespace operators {

template<class T> struct algebraic_traits {static auto one() { return T{1.0}; }};

template<class T> struct algebraic_traits<std  ::complex<T>> {static auto one() {return std  ::complex<T>{T{1}, T{0}};}};
// template<class T> struct algebraic_traits<multi::complex<T>> {static auto one() {return multi::complex<T>{T{1}, T{0}};}};

template<class X1D, class Y1D> auto operator+=(X1D&& x, Y1D const& other) BOOST_MULTI_DECLRETURN(axpy(+algebraic_traits<typename Y1D::value_type>::one(), other, std::forward<X1D>(x)))  // NOLINT(fuchsia-default-arguments-calls,readability-identifier-length) conventional name in BLAS
template<class X1D, class Y1D> auto operator-=(X1D&& x, Y1D const& other) BOOST_MULTI_DECLRETURN(axpy(-algebraic_traits<typename Y1D::value_type>::one(), other, std::forward<X1D>(x)))  // NOLINT(fuchsia-default-arguments-calls,readability-identifier-length) conventional name in BLAS

template<class X,
	std::enable_if_t<X::dimensionality == 1, int> =0>  // NOLINT(modernize-use-constraints) for C++20
auto operator*(typename X::element_type a, X const& x) {return scaled{a, x};}  // NOLINT(readability-identifier-length) conventional BLAS naming

template<class X1D, class Y1D> auto operator+(X1D const& x, Y1D const& y) -> std::decay_t<decltype(x.decay())> {auto X = x.decay(); X += y; return X;}  // NOLINT(readability-identifier-length) conventional name in BLAS
template<class X1D, class Y1D> auto operator-(X1D const& x, Y1D const& y) -> std::decay_t<decltype(x.decay())> {auto X = x.decay(); X -= y; return X;}  // NOLINT(readability-identifier-length) conventional name in BLAS

} // end namespace operators

} // end namespace boost::multi::blas

#undef BOOST_MULTI_DECLRETURN
#undef BOOST_MULTI_JUSTRETURN

#endif  // BOOST_MULTI_ADAPTORS_BLAS_AXPY_HPP

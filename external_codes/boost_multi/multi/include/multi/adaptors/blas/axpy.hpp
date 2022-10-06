// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2021 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_AXPY_HPP
#define MULTI_ADAPTORS_BLAS_AXPY_HPP

#include "../../array_ref.hpp"

#include "../../adaptors/blas/core.hpp"

#include "../../config/NODISCARD.hpp"

namespace boost::multi::blas {

using core::axpy;

template<class It1, class Size, class OutIt>
auto axpy_n(typename It1::value_type alpha, It1 first, Size n, OutIt d_first)
->decltype(axpy(n, &alpha, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n) {
	return axpy(n, &alpha, base(first) , stride(first) , base(d_first) , stride(d_first) ), d_first + n; }

template<class Context, class It1, class Size, class OutIt, class=std::enable_if_t<is_context<Context>{}>>
auto axpy_n(Context&& ctxt, typename It1::value_type alpha, It1 first, Size n, OutIt d_first)
->decltype(std::forward<Context>(ctxt).axpy(n, &alpha, first.base(), first.stride(), d_first.base(), d_first.stride()), d_first + n) {
	return std::forward<Context>(ctxt).axpy(n, &alpha, base(first) , stride(first) , base(d_first) , stride(d_first)) , d_first + n; }

template<class X1D, class Y1D, typename = decltype( std::declval<Y1D&&>()[0] = 0. )>
auto axpy(typename X1D::element alpha, X1D const& x, Y1D&& y)  // NOLINT(readability-identifier-length) conventional BLAS names
->decltype(/*axpy_n(alpha, x.begin(), x.size(), y.begin()),*/ axpy_n(alpha, x.begin(), size(x), y.begin()), std::forward<Y1D>(y)) {
	assert(size(x)==size(y)); // intel doesn't like ADL in deduced/sfinaed return types // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : bug in clang-tidy https://reviews.llvm.org/D31130
	return axpy_n(alpha, begin(x), size(x), begin(y)), std::forward<Y1D>(y);
}

template<class Context, class X1D, class Y1D, typename = decltype( std::declval<Y1D&&>()[0] = 0. )>
auto axpy(Context&& ctxt, typename X1D::element alpha, X1D const& x, Y1D&& y)  // NOLINT(readability-identifier-length) conventional BLAS names
->decltype(/*axpy_n(std::forward<Context>(ctxt), alpha, x.begin( ), x.size( ), y.begin( )),*/ std::forward<Y1D>(y)) {
	assert(size(x)==size(y)); // intel doesn't like ADL in deduced/sfinaed return types // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : bug in clang-tidy https://reviews.llvm.org/D31130
	return axpy_n(std::forward<Context>(ctxt), alpha,   begin(x),   size(x),   begin(y)), std::forward<Y1D>(y); }

template<class X1D, class Y1D>
auto axpy(X1D const& x, Y1D&& y) -> Y1D&& {  // NOLINT(readability-identifier-length) conventional BLAS names
	return axpy(+1., x, std::forward<Y1D>(y));
}

template<class Context, class X1D, class Y1D, std::enable_if_t<is_context<Context>{}> >
auto axpy(Context&& ctxt, X1D const& x, Y1D&& y) -> Y1D&& {  // NOLINT(readability-identifier-length) conventional BLAS names
	return axpy(std::forward<Context>(ctxt), +1., x, std::forward<Y1D>(y));
}

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

	template<class Other>
	friend auto operator+=(Other&& other, axpy_range const& self) -> Other&& {
		assert(other.size() == self.count_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : bug in clang-tidy https://reviews.llvm.org/D31130
		blas::axpy_n(std::forward<Context>(self.ctxt_), +self.alpha_, self.x_begin_, self.count_, other.begin());
		return std::forward<Other>(other);
	}
	template<class Other>
	friend auto operator-=(Other&& other, axpy_range const& self) -> Other&& {
		assert(other.size() == self.count_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay) : bug in clang-tidy https://reviews.llvm.org/D31130
		blas::axpy_n(std::forward<Context>(self.ctxt_), -self.alpha_, self.x_begin_, self.count_, other.begin());
		return std::forward<Other>(other);
	}
	auto operator*=(Scale s) & -> axpy_range& {alpha_ *= s; return *this;}  // NOLINT(readability-identifier-length) conventional BLAS naming
};

template<class Context, class Scalar, class X1D, class=std::enable_if_t<is_context<Context>{}>>
auto axpy(Context&& ctxt, Scalar a, X1D const& x)  // NOLINT(readability-identifier-length) conventional BLAS naming
-> axpy_range<Context, Scalar, typename X1D::const_iterator> {  // NOLINT(readability-identifier-length) conventional BLAS naming
	return {std::forward<Context>(ctxt), a, begin(x), end(x)};
}

template<class Scalar, class X1D>
auto axpy(Scalar a, X1D const& x)  // NOLINT(readability-identifier-length) conventional BLAS naming
-> axpy_range<blas::context const&, Scalar, typename X1D::const_iterator> {
	static blas::context ctxt{};
	return {ctxt, a, begin(x), end(x)};  // TODO(correaa) fix temporary
}

namespace operators {

template<class X1D, class Y1D> auto operator+=(X1D&& x, Y1D const& other) DECLRETURN(axpy(+1., other, std::forward<X1D>(x)))  // NOLINT(readability-identifier-length) conventional name in BLAS
template<class X1D, class Y1D> auto operator-=(X1D&& x, Y1D const& other) DECLRETURN(axpy(-1., other, std::forward<X1D>(x)))  // NOLINT(readability-identifier-length) conventional name in BLAS

template<class X1D, class Y1D> auto operator+(X1D const& x, Y1D const& y) -> std::decay_t<decltype(x.decay())> {auto X = x.decay(); X += y; return X;}  // NOLINT(readability-identifier-length) conventional name in BLAS
template<class X1D, class Y1D> auto operator-(X1D const& x, Y1D const& y) -> std::decay_t<decltype(x.decay())> {auto X = x.decay(); X -= y; return X;}  // NOLINT(readability-identifier-length) conventional name in BLAS

} // end namespace operators

} // end namespace boost::multi::blas
#endif

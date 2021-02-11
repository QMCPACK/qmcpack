// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_AXPY_HPP
#define MULTI_ADAPTORS_BLAS_AXPY_HPP

#include "../../adaptors/blas/core.hpp"
#include "../../config/NODISCARD.hpp"
#include "../../array_ref.hpp"

namespace boost{
namespace multi{namespace blas{

using core::axpy;

template<class T, class It1, class Size, class OutIt>
auto axpy_n(T alpha, It1 first, Size n, OutIt d_first)
->decltype(axpy(n, alpha, base(first), stride(first), base(d_first), stride(d_first)), d_first + n){
	return axpy(n, alpha, base(first), stride(first), base(d_first), stride(d_first)), d_first + n;}

template<class Context, class T, class It1, class Size, class OutIt, class=std::enable_if_t<is_context<Context>{}>>
auto axpy_n(Context&& ctxt, T alpha, It1 first, Size n, OutIt d_first)
->decltype(std::forward<Context>(ctxt).axpy(n, alpha, base(first), stride(first), base(d_first), stride(d_first)), d_first + n){
	return std::forward<Context>(ctxt).axpy(n, alpha, base(first), stride(first), base(d_first), stride(d_first)), d_first + n;}

template<class X1D, class Y1D, typename = decltype( std::declval<Y1D&&>()[0] = 0. )>
auto axpy(typename X1D::element alpha, X1D const& x, Y1D&& y)
->decltype(axpy_n(alpha, x.begin(), x.size(), y.begin()), std::forward<Y1D>(y)){assert(size(x)==size(y)); // intel doesn't like ADL in deduced/sfinaed return types
	return axpy_n(alpha, begin(x), size(x), begin(y)), std::forward<Y1D>(y);
}

template<class Context, class X1D, class Y1D, typename = decltype( std::declval<Y1D&&>()[0] = 0. )>
auto axpy(Context&& ctxt, typename X1D::element alpha, X1D const& x, Y1D&& y)
->decltype(axpy_n(std::forward<Context>(ctxt), alpha, x.begin( ), x.size( ), y.begin( )), std::forward<Y1D>(y)){assert(size(x)==size(y)); // intel doesn't like ADL in deduced/sfinaed return types
	return axpy_n(std::forward<Context>(ctxt), alpha,   begin(x),   size(x),   begin(y)), std::forward<Y1D>(y);
}

template<class X1D, class Y1D>
Y1D&& axpy(X1D const& x, Y1D&& y){return axpy(+1., x, std::forward<Y1D>(y));}

template<class Context, class X1D, class Y1D, std::enable_if_t<is_context<Context>{}> >
Y1D&& axpy(Context&& ctxt, X1D const& x, Y1D&& y){return axpy(std::forward<Context>(ctxt), +1., x, std::forward<Y1D>(y));}

template<class Context, class Scale, class ItX>
class axpy_range{
	Context ctxt_;
	Scale alpha_;
	ItX x_begin_;
	size_type count_;
public:
	axpy_range(axpy_range const&) = delete;
	axpy_range(Context ctxt, Scale alpha, ItX x_first, ItX x_last)
		: ctxt_{ctxt}, alpha_{alpha}, x_begin_{x_first}, count_{x_last - x_first}{}
	template<class Other>
	friend Other&& operator+=(Other&& other, axpy_range const& self){
		assert(other.size() == self.count_);
		blas::axpy_n(std::forward<Context>(self.ctxt_), +self.alpha_, self.x_begin_, self.count_, other.begin());
		return std::forward<Other>(other);
	}
	template<class Other>
	friend Other&& operator-=(Other&& other, axpy_range const& self){
		assert(other.size() == self.count_);
		blas::axpy_n(std::forward<Context>(self.ctxt_), -self.alpha_, self.x_begin_, self.count_, other.begin());
		return std::forward<Other>(other);
	}
	axpy_range& operator*=(Scale s)&{alpha_ *= s;}
};

template<class Context, class Scale, class X, class=std::enable_if_t<is_context<Context>{}>>
axpy_range<Context, Scale, typename X::const_iterator> axpy(Context&& ctxt, Scale a, X const& x){
	return {std::forward<Context>(ctxt), a, begin(x), end(x)};}

template<class Scale, class X>
axpy_range<blas::context const&, Scale, typename X::const_iterator> axpy(Scale a, X const& x){return {blas::context{}, a, begin(x), end(x)};}

namespace operators{

template<class X1D, class Y1D> auto operator+=(X1D&& x, Y1D const& other) DECLRETURN(axpy(+1., other, std::forward<X1D>(x)))
template<class X1D, class Y1D> auto operator-=(X1D&& x, Y1D const& other) DECLRETURN(axpy(-1., other, std::forward<X1D>(x)))

template<class X1D, class Y1D> auto operator+(X1D const& x, Y1D const& y)->std::decay_t<decltype(x.decay())>{auto X=x.decay(); X+=y; return X;}
template<class X1D, class Y1D> auto operator-(X1D const& x, Y1D const& y)->std::decay_t<decltype(x.decay())>{auto X=x.decay(); X-=y; return X;}

}


}}

}
#endif


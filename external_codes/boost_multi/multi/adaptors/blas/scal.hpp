#ifndef MULTI_ADAPTORS_BLAS_SCAL_HPP // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
#define MULTI_ADAPTORS_BLAS_SCAL_HPP
// © Alfredo A. Correa 2019-2021

#include "../blas/core.hpp"

namespace boost{
namespace multi::blas{

using core::scal;

template<class A, class It, class Size>
auto scal_n(A const& a, It first, Size count)
->decltype(scal(count, &a, first.base(), first.stride()), void()){
	       scal(count, &a, first.base(), first.stride());        }

template<class A, class It1D>
auto scal(A const& a, It1D first, It1D last)
->decltype(blas::scal_n(a, first, last - first)){
	return blas::scal_n(a, first, last - first);}

template<class A, class X1D> // don't do this: ", typename Elem = typename X1D::element_type>"
auto scal(A const& a, X1D&& x)
->decltype(blas::scal(a, x.begin(), x.end()), std::forward<X1D>(x)){
	return blas::scal(a, x.begin(), x.end()), std::forward<X1D>(x);}

template<class A>
class scal_range{
	A alpha_;
public:
	using scalar_type = A;
	explicit scal_range(A const& alpha) : alpha_{alpha}{}
	template<class X1D>
	friend auto operator*=(X1D&& x, scal_range const& self)
	->decltype(std::forward<X1D>(scal(std::declval<scalar_type const&>(), x))){
		return std::forward<X1D>(scal(self.alpha_, x));}
};

template<class A> auto scal(A const& a){return scal_range<A>{a};}

} // end namespace multi::blas
} // end namespace boost

#endif


// Copyright 2019-2026 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_TRSV_HPP
#define MULTI_ADAPTORS_BLAS_TRSV_HPP

#include "../blas/core.hpp"

#include "../blas/operations.hpp" // uplo
#include "../blas/filling.hpp"
#include "../blas/side.hpp"

#include "../../config/NODISCARD.hpp"

namespace boost::multi::blas {

enum class diagonal : char {//typename std::underlying_type<char>::type{
	unit = 'U',
	non_unit = 'N', general = non_unit
};

using core::trsv;

template<class A, std::enable_if_t<not is_conjugated<A>{}, int> =0> 
auto trsv_base(A&& a) {return base(a);}

template<class A, std::enable_if_t<    is_conjugated<A>{}, int> =0> 
auto trsv_base(A&& a) {return underlying(base(a));}

template<class A2D, class X1D>
auto trsv(filling a_nonzero_side, diagonal a_diag, A2D const& a, X1D&& x)
->decltype(trsv(static_cast<char>(flip(a_nonzero_side)), 'N', static_cast<char>(a_diag), size(x), trsv_base(a), stride(rotated(a)), trsv_base(x), stride(x)), std::forward<X1D>(x))
{
//  if(is_conjugated(x)) trsv(a_nonzero_side, a_diag, conjugated(a), conjugated(std::forward<X1D>(x)));
	{
		auto base_a = trsv_base(a);
		auto base_x = trsv_base(x);
		if(! is_conjugated<A2D>{}) {
			if     (stride(        a )==1) {trsv(static_cast<char>(flip(a_nonzero_side)), 'N', static_cast<char>(a_diag), size(x), base_a, stride(rotated(a)), base_x, stride(x));}
			else if(stride(rotated(a))==1) {trsv(static_cast<char>(     a_nonzero_side ), 'T', static_cast<char>(a_diag), size(x), base_a, stride(        a ), base_x, stride(x));}
			else {assert(0);}
		}else{
			if     (stride(        a )==1) {assert(0);}  // TODO(correaa) fallback to trsm?
			else if(stride(rotated(a))==1) {trsv(static_cast<char>(     a_nonzero_side ), 'C', static_cast<char>(a_diag), size(x), base_a, stride(        a ), base_x, stride(x));}
			else {assert(0);}
		}
	}
	return std::forward<X1D>(x);
}

template<class A2D, class X1D>
auto trsv(filling a_nonzero_side, A2D const& a, X1D&& x)
->decltype(trsv(a_nonzero_side, diagonal::general, a, std::forward<X1D>(x))) {
	return trsv(a_nonzero_side, diagonal::general, a, std::forward<X1D>(x)); }

}  // end namespace boost::multi::blas
#endif

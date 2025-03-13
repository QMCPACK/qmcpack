// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_SYRK_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_SYRK_HPP
#pragma once

#include "../blas/core.hpp"
#include "../blas/filling.hpp"
#include "../blas/numeric.hpp"

namespace boost::multi::blas {

using core::syrk;

template<class A2D, class C2D>
auto syrk(filling c_side, typename A2D::element alpha, A2D const& a, typename A2D::element beta, C2D&& c) {  // NOLINT(readability-identifier-length) BLAS naming
	//->decltype(syrk('\0', '\0', size(c), size(a), alpha, base(a), stride(rotated(a)), beta, base(c), stride(c)), std::forward<C2D>(c)){
	assert(size(c) == size(rotated(c)));
	if(stride(a) == 1) {
		if(stride(c) == 1) {
			syrk(flip(c_side) == filling::upper ? 'L' : 'U', 'N', size(c), size(a), &alpha, base(a), stride(rotated(a)), &beta, base(c), stride(rotated(c)));
		} else {
			syrk(c_side == filling::upper ? 'L' : 'U', 'N', size(c), size(rotated(a)), &alpha, base(a), stride(rotated(a)), &beta, base(c), stride(c));
		}
	} else {
		if(stride(c) == 1) {
			syrk(flip(c_side) == filling::upper ? 'L' : 'U', 'T', size(c), size(rotated(a)), &alpha, base(a), stride(a), &beta, base(c), stride(rotated(c)));
		} else {
			syrk(c_side == filling::upper ? 'L' : 'U', 'T', size(c), size(rotated(a)), &alpha, base(a), stride(a), &beta, base(c), stride(c));
		}
	}
	return std::forward<C2D>(c);
}

template<typename AA, class A2D, class C2D>
auto syrk(filling c_side, AA alpha, A2D const& a, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
	-> decltype(syrk(c_side, alpha, a, 0.0, std::forward<C2D>(c))) {
	return syrk(c_side, alpha, a, 0.0, std::forward<C2D>(c));
}

// template<typename AA, class A2D, class C2D>
// auto syrk(AA alpha, A2D const& a, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
//  -> decltype(syrk(filling::upper, alpha, a, syrk(filling::lower, alpha, a, std::forward<C2D>(c)))) {
//  return syrk(filling::upper, alpha, a, syrk(filling::lower, alpha, a, std::forward<C2D>(c)));
// }

// template<typename AA, class A2D, class Ret = typename A2D::decay_type>
// [[nodiscard]]  // ("because input argument is const")
// // this decay in the return type is important
// auto  // NOLINTNEXTLINE(readability-identifier-length) BLAS naming
// syrk(AA alpha, A2D const& a) -> std::decay_t<decltype(syrk(alpha, a, Ret({size(a), size(a)}, get_allocator(a))))> {
//  return syrk(alpha, a, Ret({size(a), size(a)}, get_allocator(a)));
// }

// template<class A2D>
// [[nodiscard]] auto syrk(A2D const& A)  // NOLINT(readability-identifier-length) BLAS naming
//  -> decltype(syrk(1.0, A)) {
//  return syrk(1.0, A);
// }

}  // end namespace boost::multi::blas
#endif

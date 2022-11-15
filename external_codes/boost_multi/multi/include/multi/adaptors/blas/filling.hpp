// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_FILLING_HPP
#define MULTI_ADAPTORS_BLAS_FILLING_HPP

#include "multi/array_ref.hpp"

#include    "../blas/core.hpp"
#include    "../blas/operations.hpp"

namespace boost::multi::blas {

enum class filling : char {
	lower = 'U',
	upper = 'L'
};

[[maybe_unused]] static constexpr filling U = filling::upper;  // NOLINT(readability-identifier-length) BLAS naming
[[maybe_unused]] static constexpr filling L = filling::lower;  // NOLINT(readability-identifier-length) BLAS naming

inline auto flip(filling side) -> filling {
	switch(side) {
		case filling::lower: return filling::upper;
		case filling::upper: return filling::lower;
	} __builtin_unreachable();
}

inline auto operator-(filling side) -> filling {return flip(side);}
inline auto operator+(filling side) -> filling {return side;}

template<class A2D, std::enable_if_t<is_conjugated<A2D>{}, int> =0>
auto detect_triangular_aux(A2D const& array, std::false_type /*false*/) -> filling {
	{
		for(auto i = size(array); i != 0; --i) {  // NOLINT(altera-id-dependent-backward-branch)
			auto const asum_up = blas::asum(begin(array[i-1])+i, end(array[i-1]));
			if(std::isnan(asum_up)) {return filling::lower;}
			if(asum_up !=0.       ) {return filling::upper;}

			auto const asum_lo = blas::asum(begin(rotated(array)[i-1])+i, end(rotated(array)[i-1]));
			if(std::isnan(asum_lo)) {return filling::upper;}
			if(asum_lo != 0.      ) {return filling::lower;}
		}
	}
	return filling::lower;
}

template<class A2D>
auto detect_triangular(A2D const& /*array*/) -> filling;

template<class A2D, std::enable_if_t<is_conjugated<A2D>{}, int> =0>
auto detect_triangular_aux(A2D const& array) -> filling {
	return flip(detect_triangular(hermitized(array)));
}

template<class A2D>
auto detect_triangular(A2D const& array) -> filling {
#if defined(__cpp_if_constexpr)
	if constexpr(not is_conjugated<A2D>{}) {
		using blas::asum;
		for(auto i = size(array); i != 0; --i) {
			auto const asum_up = asum(A[i-1]({i, array[i-1].size()}));
			if(std::isnan(asum_up)) {return filling::lower;}
			if(asum_up!=0.        ) {return filling::upper;}

			auto const asum_lo = asum(rotated(array)[i-1]({i, rotated(array)[i-1].size()}));
			if(std::isnan(asum_lo)) {return filling::upper;}
			if(asum_lo != 0.      ) {return filling::lower;}
		}
		return filling::lower;
	} else {
		return flip(detect_triangular(hermitized(array)));
	}
#else
	return detect_triangular_aux(A);//, is_conjugated<A2D>{});//std::integral_constant<bool, not is_hermitized<A2D>()>{});
#endif
}

}  // end namespace boost::multi::blas

#endif

// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_FILLING_HPP
#define MULTI_ADAPTORS_BLAS_FILLING_HPP

#include <multi/array_ref.hpp>

#include <multi/adaptors/blas/core.hpp>
#include <multi/adaptors/blas/operations.hpp>

namespace boost::multi::blas {

enum class filling : char {
	lower = 'U',
	upper = 'L'
};

inline auto flip(filling side) -> filling {
	switch(side) {
		case filling::lower: return filling::upper;
		case filling::upper: return filling::lower;
	} __builtin_unreachable();  // LCOV_EXCL_LINE
}

inline auto operator-(filling side) -> filling {return flip(side);}
inline auto operator+(filling side) -> filling {return side;}

}  // end namespace boost::multi::blas

#endif

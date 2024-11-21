// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_FILLING_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_FILLING_HPP

#include <boost/multi/array_ref.hpp>

#include <boost/multi/adaptors/blas/core.hpp>
#include <boost/multi/adaptors/blas/operations.hpp>

namespace boost::multi::blas {

enum class filling : char {
	lower = 'U',
	upper = 'L'
};

inline auto flip(filling side) -> filling {
	switch(side) {
		case filling::lower: return filling::upper;
		case filling::upper: return filling::lower;
	}  // __builtin_unreachable();  // LCOV_EXCL_LINE
	return {};
}

inline auto operator-(filling side) -> filling {return flip(side);}
inline auto operator+(filling side) -> filling {return side;}

}  // end namespace boost::multi::blas

#endif

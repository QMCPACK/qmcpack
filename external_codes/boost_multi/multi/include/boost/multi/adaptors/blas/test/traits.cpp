// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/blas/traits.hpp>
#include <boost/multi/array.hpp>

#include <complex>

namespace multi = boost::multi;
namespace blas  = multi::blas;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_traits_simple_array) {
		multi::array<double, 2> const arr;
		BOOST_TEST( arr.empty() );
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_traits) {
		static_assert(blas::is_d<double>{});
		static_assert(blas::is_s<float>{});

		static_assert(blas::is_c<std::complex<float>>{});
		static_assert(blas::is_z<std::complex<double>>{});
	}
	return boost::report_errors();
}

// Copyright 2018-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array

#include <boost/core/lightweight_test.hpp>

#include <complex>  // for complex
#include <string>   // for operator==, string

namespace multi = boost::multi;

namespace {
inline auto what_is(multi::array<double, 2> const& /*arr*/) { return std::string{"real"}; }                   // NOLINT(fuchsia-default-arguments-calls)
inline auto what_is(multi::array<std::complex<double>, 2> const& /*arr*/) { return std::string{"complex"}; }  // NOLINT(fuchsia-default-arguments-calls)
}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(multi_array_overload)
	{
		multi::array<double, 2> const               real_A({10, 20});
		multi::array<std::complex<double>, 2> const cplx_A({10, 20});

		std::string const real_str    = what_is(real_A);
		std::string const complex_str = what_is(cplx_A);

		BOOST_TEST( real_str    == "real"    );
		BOOST_TEST( complex_str == "complex" );
	}
	return boost::report_errors();
}

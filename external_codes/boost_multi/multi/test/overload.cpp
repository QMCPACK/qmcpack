// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #if defined(__clang__)
//  #pragma clang diagnostic push
//  #pragma clang diagnostic ignored "-Wunknown-warning-option"
//  #pragma clang diagnostic ignored "-Wconversion"
//  #pragma clang diagnostic ignored "-Wextra-semi-stmt"
//  #pragma clang diagnostic ignored "-Wold-style-cast"
//  #pragma clang diagnostic ignored "-Wsign-conversion"
//  #pragma clang diagnostic ignored "-Wswitch-default"
//  #pragma clang diagnostic ignored "-Wundef"
// #elif defined(__GNUC__)
//  #pragma GCC diagnostic push
//  #if (__GNUC__ > 7)
//      #pragma GCC diagnostic ignored "-Wcast-function-type"
//  #endif
//  #pragma GCC diagnostic ignored "-Wconversion"
//  #pragma GCC diagnostic ignored "-Wold-style-cast"
//  #pragma GCC diagnostic ignored "-Wsign-conversion"
//  #pragma GCC diagnostic ignored "-Wundef"
// #endif

// #ifndef BOOST_TEST_MODULE
//  #define BOOST_TEST_MAIN
// #endif

// #include <boost/test/included/unit_test.hpp>

// #if defined(__clang__)
//  #pragma clang diagnostic pop
// #elif defined(__GNUC__)
//  #pragma GCC diagnostic pop
// #endif

#include <boost/multi/array.hpp>  // for array

#include <complex>  // for complex
#include <string>   // for operator==, string

namespace multi = boost::multi;

inline auto what_is(multi::array<double, 2> const& /*arr*/) { return std::string{"real"}; }                   // NOLINT(fuchsia-default-arguments-calls)
inline auto what_is(multi::array<std::complex<double>, 2> const& /*arr*/) { return std::string{"complex"}; }  // NOLINT(fuchsia-default-arguments-calls)

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(multi_array_overload) {
	multi::array<double, 2> const               real_A({10, 20});
	multi::array<std::complex<double>, 2> const cplx_A({10, 20});

	std::string const real_str    = what_is(real_A);
	std::string const complex_str = what_is(cplx_A);

	BOOST_TEST( real_str    == "real"    );
	BOOST_TEST( complex_str == "complex" );
}
return boost::report_errors();}

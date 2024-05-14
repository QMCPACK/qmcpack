// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi overload resolution"  // title NOLINT(cppcoreguidelines-macro-usage)
#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include <complex>

namespace multi = boost::multi;

inline auto what_is(multi::array<double              , 2> const& /*arr*/) { return std::string{"real"}; }  // std::string NOLINT(fuchsia-default-arguments-calls)
inline auto what_is(multi::array<std::complex<double>, 2> const& /*arr*/) { return std::string{"complex"}; }  // std::string NOLINT(fuchsia-default-arguments-calls)

BOOST_AUTO_TEST_CASE(multi_array_overload) {
	multi::array<double, 2> const               real_A({10, 20});
	multi::array<std::complex<double>, 2> const cplx_A({10, 20});

	std::string const real_str    = what_is(real_A);
	std::string const complex_str = what_is(cplx_A);

	BOOST_REQUIRE( real_str    == "real"    );
	BOOST_REQUIRE( complex_str == "complex" );
}

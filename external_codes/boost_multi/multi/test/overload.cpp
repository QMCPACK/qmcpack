// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi overload resolution"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>

namespace multi = boost::multi;

inline auto what_is(multi::array<             double , 2> const& /*arr*/){return std::string{"real"}   ;}
inline auto what_is(multi::array<std::complex<double>, 2> const& /*arr*/){return std::string{"complex"};}

BOOST_AUTO_TEST_CASE(multi_array_range_section) {
	multi::array<             double , 2> real_A({10, 20});
	multi::array<std::complex<double>, 2> cplx_A({10, 20});

	std::string real_str    = what_is(real_A);
	std::string complex_str = what_is(cplx_A);

	BOOST_REQUIRE( real_str    == "real"    );
	BOOST_REQUIRE( complex_str == "complex" );
}


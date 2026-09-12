#define BOOST_TEST_MODULE "C++ Unit Tests for Multi TotalView adaptor"
#define BOOST_TEST_DYN_LINK

// #include <boost/test/unit_test.hpp>

#include "multi/array.hpp"
#include "multi/utility.hpp"

#include "../../../adaptors/totalview.hpp"

#include <algorithm>  // transform
#include <complex>
#include <iostream>
#include <numeric>  // iota

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_1d) {

	std::vector<int> V = {10, 20, 30};

	multi::array<double, 1> const A     = {1.0, 2.0, 3.0, 4.0, 5.0};
	auto&&                        Apart = A({1, 3});

	multi::array<double, 2> const B = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
	};

	double sum = 0.0;
	for(auto i : A.extension()) {
		sum += A[i];
	}

	BOOST_REQUIRE( sum == 15.0 );
	BOOST_REQUIRE( B[1][0] == 4.0 );
}

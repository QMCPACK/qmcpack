// Copyright 2023-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

namespace multi = boost::multi;

auto main() -> int {
	multi::array<int, 2> arr = {
		{ 0,  1,  2,  3,  4},
		{ 5,  6,  7,  8,  9},
		{10, 11, 12, 13, 14},
		{15, 16, 17, 18, 19},
	};

	if(arr[2][3] != 13) {
		return 1;
	}
	return 0;
}

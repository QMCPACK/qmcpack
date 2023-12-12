// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#include <multi/array.hpp>

namespace multi = boost::multi;

auto main() -> int {
	multi::array<double, 2> arr = {
		{ 0.0,  1.0,  2.0,  3.0,  4.0},
		{ 5.0,  6.0,  7.0,  8.0,  9.0},
		{10.0, 11.0, 12.0, 13.0, 14.0},
		{15.0, 16.0, 17.0, 18.0, 19.0},
	};

	if(arr[2][3] != 13.0) {
		return 1;
	}
	return 0;
}

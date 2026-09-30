// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, dynamic_array, num_elements

#include <boost/core/lightweight_test.hpp>

#include <utility>  // for std::move

namespace multi = boost::multi;

namespace {
auto fun(multi::subarray<int, 2, boost::multi::detail::offset_ptr<int>>& arr2d) -> int {
	return arr2d[1][1];
}

auto gun(multi::subarray<int, 2> arr2d) -> int {
	arr2d[1][0] = 0;
	return arr2d[1][1];
}
}  // namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	multi::inplace_array<int[4][4]> a2d = {
		{1, 2},
		{3, 4}
	};

	BOOST_TEST( a2d.size() == 2 );
	BOOST_TEST( a2d[1][1] == 4 );

	// UB probe: form a pointer that is *past* one-past-end of a transposed
	{
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		multi::inplace_array<int[2][2]> b2d = {
			{1, 2},
			{3, 4}
		};
		auto const end = (~b2d)[1].end();    // UB: forms a pointer past one-past-end
		BOOST_TEST( end != (~b2d)[1].begin() );  // observe the pointer so the optimizer can't drop the arithmetic
	}

	fun(a2d);

	multi::subarray<int, 2, boost::multi::detail::offset_ptr<int>>& a2ds = a2d;

	gun(std::move(a2ds));

	return boost::report_errors();
}

// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <cstddef>   // for std::ptrdiff_t
#include <iterator>  // for std::distance
#include <numeric>   // for std::iota

namespace multi = boost::multi;

namespace {

// Case 1: Strided 1D view where stride does not divide extent.
void test_strided_1d_index_access() {
	multi::array<int, 1> arr(9);
	std::iota(arr.begin(), arr.end(), 0);

	auto strided = arr.strided(3);
	BOOST_TEST( strided.size() == 3 );

	// Index-based access works correctly:
	BOOST_TEST( strided[0] == 0 );
	BOOST_TEST( strided[1] == 3 );
	BOOST_TEST( strided[2] == 6 );

	// But iteration via begin()/end() would be an infinite loop,
}

// Case 2: Transposed 2D array — subarray .end() overshoots allocation.
void test_transposed_subarray_overshoot() {
	multi::array<int, 2> arr({2, 5});
	std::iota(arr.elements().begin(), arr.elements().end(), 0);

	auto transposed = arr.transposed();  // 5x2 view
	BOOST_TEST( transposed.size() == 5 );

	// Row 4 of the transposed view: elements at positions 4 and 9
	auto row4 = transposed[4];
	BOOST_TEST( row4[0] == 4 );
	BOOST_TEST( row4[1] == 9 );
	BOOST_TEST( row4.size() == 2 );

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif
	// row4.end() forms a pointer at base + 4 + 5*2 = base + 14,
	// but allocation ends at base + 10.  This is 4 elements past the allocation.
	// UB under [expr.add], but works on all tested platforms.
	// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	std::ptrdiff_t const overshoot = (row4.base() + (row4.stride() * row4.size())) - (arr.base() + arr.num_elements());

	BOOST_TEST( overshoot == 4 );  // 14 - 10 = 4
#ifdef __clang__
#pragma clang diagnostic pop
#endif

	// Iteration works because stride divides the subarray's nelems:
	auto end = row4.end();
	BOOST_TEST( std::distance(row4.begin(), end) == 2 );
}

// Case 3: Strided 1D view where stride divides extent — no problem.
void test_strided_1d_divisible() {
	multi::array<int, 1> arr(12);
	std::iota(arr.begin(), arr.end(), 0);

	auto strided = arr.strided(3);
	BOOST_TEST( strided.size() == 4 );

	// This is fine: end pointer = base + 12 = one-past-end of allocation.
	int sum = 0;
	for(auto it = strided.begin(); it != strided.end(); ++it) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch,modernize-loop-convert)
		sum += *it;
	}
	BOOST_TEST( sum == 0 + 3 + 6 + 9 );
}

}  // end namespace

auto main() -> int {  // NOLINT(bugprone-exception-escape)
	test_strided_1d_index_access();
	test_transposed_subarray_overshoot();
	test_strided_1d_divisible();

	{
		multi::array<double, 2> arr = {
			{1.0, 2.0, 3.0},
			{3.0, 4.0, 5.0},
			{6.0, 7.0, 8.0},
			{9.0, 9.0, 9.0}
		};

		BOOST_TEST( arr.strided(2).size() == 2 );
	}
	// {
	// 	multi::array<double, 2> arr = {
	// 		{1.0, 2.0, 3.0},
	// 		{3.0, 4.0, 5.0},
	// 		{6.0, 7.0, 8.0},
	// 	};

	// 	BOOST_TEST( arr.strided(2).size() == 2 );
	// }

	return boost::report_errors();
}

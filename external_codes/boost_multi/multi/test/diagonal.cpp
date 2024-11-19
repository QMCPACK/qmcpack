// Copyright 2023-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, layout_t, static_array

#include <algorithm>   // for transform
#include <functional>  // for plus  // IWYU pragma: keep
#include <numeric>     // for accumulate

namespace multi = boost::multi;

template<class Array2D>
auto trace_with_indices(Array2D const& arr) {
	typename Array2D::element_type sum{0};
	for(auto i : extension(arr)) {  // NOLINT(altera-unroll-loops) testing loops
		sum += arr[i][i];
	}
	return sum;
}

template<class Array2D>
auto trace_with_diagonal(Array2D const& arr) {
	typename Array2D::element_type sum{0};
	for(auto aii : arr.diagonal()) {  // NOLINT(altera-unroll-loops) testing loops
		sum += aii;
	}
	return sum;
}

template<class Array2D>
auto trace_with_accumulate(Array2D const& arr) {
	return std::accumulate(arr.diagonal().begin(), arr.diagonal().end(), static_cast<typename Array2D::element_type>(0));
}

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(trace_test) {
	using int_element = multi::index;
	multi::array<int_element, 2> arr({5, 5}, 0);

	auto [is, js] = extensions(arr);

	// NOLINTNEXTLINE(altera-unroll-loops) testing loops
	for(auto i : is) {
		for(auto j : js) {  // NOLINT(altera-unroll-loops) testing loops
			arr[i][j] = 10 * i + j;
		}
	}

	auto tr = trace_with_diagonal(arr);

	BOOST_TEST( tr == 00 + 11 + 22 + 33 + 44 );

	BOOST_TEST( trace_with_diagonal(arr) == trace_with_indices(arr) );
	BOOST_TEST( trace_with_diagonal(arr) == trace_with_accumulate(arr) );
}

BOOST_AUTO_TEST_CASE(broadcasted) {
	multi::array<int, 2> const arr = {
		{0, 1,  2},
		{4, 5,  6},
		{8, 9, 10},
	};

	BOOST_TEST( arr.diagonal().begin() != arr.diagonal().end() );
	BOOST_TEST( arr.diagonal().end() - arr.diagonal().begin() == 3 );
}

BOOST_AUTO_TEST_CASE(broadcasted) {
	multi::array<int, 2> const arr = {
		{0, 1,  2,  3},
		{4, 5,  6,  7},
		{8, 9, 10, 11},
	};

	auto const& a3D = arr.broadcasted();

	BOOST_TEST( &a3D[0][2][1] == &arr[2][1] );
	BOOST_TEST( &a3D[1][2][1] == &arr[2][1] );

	{
		auto const& arr_instance = a3D[0];
		BOOST_TEST( &arr_instance[3][1] == &arr[3][1] );
	}
	{
		auto const& arr_instance = a3D[99];
		BOOST_TEST( &arr_instance[3][1] == &arr[3][1] );
	}
	{
		auto const& arr_instance = a3D[-99];
		BOOST_TEST( &arr_instance[3][1] == &arr[3][1] );
	}
	{
		auto const& a3D_self = a3D();
		BOOST_TEST( &a3D_self[ 4][3][1] == &arr[3][1] );
		BOOST_TEST( &a3D_self[99][3][1] == &arr[3][1] );
	}
}

BOOST_AUTO_TEST_CASE(broadcast_1D) {
	multi::array<int, 1> const arr = {0, 1, 2, 3};

	auto const& a2D = arr.broadcasted();

	BOOST_TEST( &a2D[0][2] == &arr[2] );
	BOOST_TEST( &a2D[1][2] == &arr[2] );
}

BOOST_AUTO_TEST_CASE(broadcast_0D) {
	multi::array<int, 1>       arr = {0, 1, 2, 3};
	multi::array<int, 0> const vv(2);

	auto const& v1D = vv.broadcasted();

	BOOST_TEST( &v1D[0] == vv.base() );
	BOOST_TEST( &v1D[1] == vv.base() );

	multi::array<int, 1> r1D({4}, 0);
	std::transform(arr.begin(), arr.end(), v1D.begin(), r1D.begin(), std::plus<>{});

	BOOST_TEST( r1D[3] == arr[3] + 2 );

	std::transform(arr.begin(), arr.end(), v1D.begin(), arr.begin(), [](auto, auto ve) { return ve; });
	BOOST_TEST( arr[3] == 2 );
}

return boost::report_errors();}

// Copyright 2023-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <numeric>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

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

BOOST_AUTO_TEST_CASE(trace_test) {
	using int_element = multi::index;
	multi::array<int_element, 2> arr({5, 5}, 0);

	auto [is, js] = extensions(arr);
	for(auto i : is) {  // NOLINT(altera-unroll-loops) testing loops
		for(auto j : js) {  // NOLINT(altera-unroll-loops) testing loops
			arr[i][j] = 10 * i + j;
		}
	}

	auto tr = trace_with_diagonal(arr);

	BOOST_REQUIRE( tr == 00 + 11 + 22 + 33 + 44 );

	BOOST_REQUIRE( trace_with_diagonal(arr) == trace_with_indices(arr) );
	BOOST_REQUIRE( trace_with_diagonal(arr) == trace_with_accumulate(arr) );
//  BOOST_REQUIRE( trace_with_diagonal(arr) == trace_with_reduce(arr) );
}

BOOST_AUTO_TEST_CASE(broadcasted) {
	multi::array<int, 2> const arr = {
		{0, 1, 2, 3},
		{4, 5, 6, 7},
		{8, 9, 10, 11},
	};

	auto const& a3D = arr.broadcasted();

	BOOST_TEST( &a3D[0][2][1] == &arr[2][1] );
	BOOST_TEST( &a3D[1][2][1] == &arr[2][1] );

	{
		auto const& arr_instance = a3D[0];
		BOOST_REQUIRE( &arr_instance[3][1] == &arr[3][1] );
	}
	{
		auto const& arr_instance = a3D[99];
		BOOST_REQUIRE( &arr_instance[3][1] == &arr[3][1] );
	}
	{
		auto const& arr_instance = a3D[-99];
		BOOST_REQUIRE( &arr_instance[3][1] == &arr[3][1] );
	}
	{
		auto const& a3D_self = a3D();
		BOOST_TEST( &a3D_self[ 4][3][1] == &arr[3][1] );
		BOOST_TEST( &a3D_self[99][3][1] == &arr[3][1] );
	}
	{
		// [[maybe_unused]] auto const& a3D_finite = a3D({0, 9});
		// BOOST_TEST( &a3D_finite[ 4][3][1] == &arr[3][1] );
		// BOOST_TEST( &a3D_finite[99][3][1] == &arr[3][1] );
	}

//  BOOST_REQUIRE( a3D_finite.size() == 5 );
//  BOOST_REQUIRE( a3D_finite.begin() + 5 == a3D_finite.end() );
}

BOOST_AUTO_TEST_CASE(broadcast_1D) {
	multi::array<int, 1> const arr = {0, 1, 2, 3};

	auto const& a2D = arr.broadcasted();

	BOOST_TEST( &a2D[0][2] == &arr[2] );
	BOOST_TEST( &a2D[1][2] == &arr[2] );
}

BOOST_AUTO_TEST_CASE(broadcast_0D) {
	multi::array<int, 1> arr = {0, 1, 2, 3};
	multi::array<int, 0> const vv(2);

	auto const& v1D = vv.broadcasted();

	BOOST_TEST( &v1D[0] == vv.base() );
	BOOST_TEST( &v1D[1] == vv.base() );

	multi::array<int, 1> r1D({4}, 0);
	std::transform(arr.begin(), arr.end(), v1D.begin(), r1D.begin(), std::plus<>{});

	BOOST_TEST( r1D[3] == arr[3] + 2 );

	std::transform(arr.begin(), arr.end(), v1D.begin(), arr.begin(), [](auto, auto ve) {return ve;});
	BOOST_TEST( arr[3] == 2 );

	// std::copy_n(v1D.begin(), arr.size(), arr.begin());
	// BOOST_TEST( arr[3] == 2 );
}

// Copyright 2019-2023 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <array>
#include <complex>

// Suppress warnings from boost.test
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wundef"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d) {
	{
		std::vector<double> const vec = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_REQUIRE( vec[1] == 2. );
	}
	{
		multi::static_array<int, 1> arr = {12, 34, 56};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 56 );
	}
	{
		multi::static_array<int, 1> const arr = {12, 34, 56};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 56 );
	}
	{
		auto const il = {1.2, 3.4, 5.6};

		multi::static_array<double, 1> const arr(il);
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == il.begin()[2] );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
	{
		auto const il = {1.2, 3.4, 5.6};

		multi::static_array<double, 1> const arr(begin(il), end(il));
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == il.begin()[2] );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
	{
		multi::static_array<int, 1> const arr = {12, 34, 56};
		BOOST_TEST_REQUIRE( size(arr) == 3 );
		BOOST_TEST_REQUIRE( arr[2] == 56 );
		BOOST_TEST_REQUIRE(( arr == multi::static_array<int, 1>{12, 34, 56} ));
		BOOST_TEST_REQUIRE(( arr == decltype(arr){12, 34, 56} ));
	}
	{
		auto const values = {12, 34, 56};

		multi::array<int, 1> const arr(values.begin(), values.end());
		BOOST_TEST_REQUIRE( size(arr) == 3 );
		BOOST_TEST_REQUIRE( arr[2] == 56 );
	}
	{
		multi::array<int, 1> const arr = {12, 34, 56};

		BOOST_TEST_REQUIRE( size(arr) == 3 );
		BOOST_TEST_REQUIRE( arr[2] == 56 );

		BOOST_TEST_REQUIRE(( arr == multi::array<int, 1>{12, 34, 56} ));
		BOOST_TEST_REQUIRE(( arr == decltype(arr){12, 34, 56} ));
		BOOST_TEST_REQUIRE(( arr == decltype(arr)::decay_type({12, 34, 56}) ));
	}
	{
		std::array<int, 3> const stdarr = {
			{11, 22, 33},
		};
		using multi::num_elements;
		BOOST_TEST_REQUIRE( num_elements(stdarr) == 3 );

		using std::begin;
		using std::end;
		multi::static_array<double, 1> const arr(begin(stdarr), end(stdarr));
		BOOST_TEST_REQUIRE( size(arr) == 3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d_ctad) {
#if defined(__cpp_deduction_guides) && !defined(__NVCC__)
#if !defined(__circle_build__) || (__circle_build__ > 200 )  // crashes circle 187-200 in docker
	{
		multi::static_array const arr = {12, 34, 56};
		BOOST_TEST_REQUIRE( size(arr) == 3 );
		BOOST_TEST_REQUIRE( arr[2] == 56 );
		BOOST_TEST_REQUIRE(( arr == multi::static_array{12, 34, 56} ));
	}
#endif
	{
		multi::array arr({12, 34, 56});
		BOOST_TEST_REQUIRE( size(arr) == 3 );
		BOOST_TEST_REQUIRE( arr[2] == 56 );
		BOOST_TEST_REQUIRE(( arr == multi::array({12, 34, 56}) ));
	}
#endif
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_array) {
#if defined(__INTEL_COMPILER) || (defined(__clang__) && (__clang_major__ >= 10))  // doesn't work on gcc
	{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-designator"
		//      double const a[] = { [8] = 8.0, 9.0, 10.0 };
		std::array<double, 11> const stdarr = {
			{[8] = 8.0, 9.0, 10.0},
		};
#pragma GCC diagnostic pop
		multi::array<double, 1> arr = stdarr;
		BOOST_REQUIRE( arr.size() == 11 );
		BOOST_REQUIRE( arr[9] == 9.0 );
	}
#endif
}

BOOST_AUTO_TEST_CASE(multi_initialize_from_carray_1d) {
	{
		multi::static_array<int, 1> const arr = {11, 22, 33};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[1] == 22 );
	}
	{
#if defined(__cpp_deduction_guides) && ! defined(__NVCC__)
//      multi::array arr = {{1.1, 2.2, 3.3}};
//      static_assert( decltype(arr)::dimensionality == 1 , "!");
//      BOOST_REQUIRE( size(arr)==3 && arr[1] == 2.2 );
#endif
	}
	{
		std::array<double, 3> stdarr = {
			{1.1, 2.2, 3.3}
		};
		multi::array<double, 1> const arr(begin(stdarr), end(stdarr));
		BOOST_REQUIRE(( arr == decltype(arr){1.1, 2.2, 3.3} ));
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_2d) {
	{
		multi::static_array<double, 2> const arr = {
			{ 1.2,  2.4, 3.6, 8.9},
			{11.2, 34.4, 5.6, 1.1},
			{15.2, 32.4, 5.6, 3.4},
		};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( size(arr[0]) == 4 );
		BOOST_REQUIRE(( arr == decltype(arr){
			{ 1.2,  2.4, 3.6, 8.9},
			{11.2, 34.4, 5.6, 1.1},
			{15.2, 32.4, 5.6, 3.4},
		}));
	}
	{
		multi::array<int, 2> const arr = {
			{ 12,  24, 36},
			{112, 344, 56},
			{152, 324, 56},
		};
		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( size(arr[0]) == 3 );
		BOOST_TEST( arr[1][1] == 344 );
	}
	{
		multi::array<int, 2> arr = {
			{ 12,  24, 36},
			{112, 344, 56},
			{152, 324, 56},
		};

		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( size(arr[0]) == 3 );
		BOOST_TEST( arr[1][1] == 344 );

		arr = {
			{100,  10,  20},
			{100, 110, 120},
			{200, 210, 220},
		};
		BOOST_TEST( arr[1][2] == 120 );
	}
	{
		multi::array<int, 1> vec;
		vec = {40, 55};
		BOOST_TEST( size(vec) == 2 );
		BOOST_TEST( vec[1] == 55 );
	}
	{
		std::array<std::array<int, 2>, 3> const nested = {
			{{{12, 24}}, {{112, 344}}, {{152, 324}}}
		};

		using std::begin;
		using std::end;

		multi::static_array<double, 2> arr(begin(nested), end(nested));

		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( size(arr[0]) == 2 );
		BOOST_TEST( arr[1][0] == 112 );
	}
	{
		std::array<std::array<int, 2>, 3> const nested = {
			{{{12, 24}}, {{112, 344}}, {{152, 324}}}
		};
		multi::static_array<int, 2> const arr(std::begin(nested), std::end(nested));

		BOOST_TEST((
			arr == multi::array<int, 2> {{
				{{ 12,  24}},
				{{112, 344}},
				{{152, 324}}
			}}
		));

		BOOST_TEST(!( arr != multi::array<int, 2>{
				{ 12,  24},
				{112, 344},
				{152, 324},
			}
		));
		BOOST_TEST((
			arr == decltype(arr){
				{ 12,  24},
				{112, 344},
				{152, 324},
			}
		));
	}
	{
		std::array<std::array<int, 2>, 3> nested = {
			{{{10, 20}},
			 {{20, 40}},
			 {{30, 60}}},
		};
		multi::array<int, 2> arr(begin(nested), end(nested));
		BOOST_TEST( num_elements(arr) == 6 );
		BOOST_TEST( arr[2][1] == 60 );
	}
	{
		using complex = std::complex<double>;

		complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imaginary unit

		multi::array<complex, 2> arr = {
			{2.0 + 1.0 * I, 1.0 + 3.0 * I, 1.0 + 7.0 * I},
			{3.0 + 4.0 * I, 4.0 + 2.0 * I, 0.0 + 0.0 * I},
		};
		BOOST_REQUIRE( arr[1][1] == 4.0 + 2.0*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_static_array_initializer_list) {
	multi::static_array<std::complex<double>, 2> SA = {
		{{1.0, 0.0}, {2.0, 0.0}},
		{{3.0, 0.0}, {4.0, 0.0}},
	};
	BOOST_REQUIRE( SA[1][1] == 4.0 );
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d) {
	multi::array<int, 3> const arr = {
		{ {12, 100},  {24, 10}},
		{ {112, 30}, {344, 40}},
		{{152, 990}, {324, 20}},
	};
	BOOST_REQUIRE( arr[1][1][0] == 344 );
	BOOST_REQUIRE( arr[1][1][1] ==  40 );
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d_string) {
	{
		using std::string;

		// NOLINTBEGIN(fuchsia-default-arguments-calls)
		multi::array<string, 3> B3 = {
			{{"000", "001", "002"}, {"010", "011", "012"}},
			{{"100", "101", "102"}, {"110", "111", "112"}},
		};
		// NOLINTEND(fuchsia-default-arguments-calls)

		BOOST_REQUIRE( num_elements(B3)==12 && B3[1][0][1] == "101" );
	}
}

#if defined(__cpp_deduction_guides) && ! defined(__NVCC__)
BOOST_AUTO_TEST_CASE(initializer_list_1d_static) {
#if !defined(__circle_build__) || (__circle_build__ > 200 )  // crashes circle 187-200 in docker
	{
		multi::static_array arr({1.0, 2.0, 3.0});
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( size(arr) == 3 && num_elements(arr) == 3 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 && num_elements(arr)==3 && arr[1] == 2.0 );
		static_assert(typename decltype(arr)::rank{} == 1);
	}
#endif
}

BOOST_AUTO_TEST_CASE(initializer_list_1d) {
	{
		multi::array arr({1.0, 2.0, 3.0});
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( size(arr) == 3 && num_elements(arr) == 3 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 && num_elements(arr)==3 && arr[1] == 2.0 );
		static_assert(typename decltype(arr)::rank{} == 1);
	}
	{
		multi::array arr({1.0, 2.0});
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( size(arr) == 2 && num_elements(arr) == 2 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 && num_elements(arr) == 2 && arr[1] == 2.0 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 );
	}
	{
		multi::array arr({0, 2});  //  multi::array arr = {0, 2}; not working with CTAD
		static_assert(std::is_same_v<decltype(arr)::element_type, int>);
		BOOST_REQUIRE( size(arr) == 2 && num_elements(arr) == 2 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 && num_elements(arr) == 2 && arr[1] == 2.0 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 );
	}
	{
		multi::array arr({9.0});  // multi::array arr = {9.0}; not working with CTAD
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 && num_elements(arr)==1 && arr[0]==9.0 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 );
	}
	{
		multi::array arr({9});  // multi::array arr = {9}; not working with CTAD
		static_assert(std::is_same_v<decltype(arr)::element_type, int>);
		BOOST_REQUIRE( size(arr) == 1 && num_elements(arr) == 1 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 );
		BOOST_REQUIRE( num_elements(arr) == 1 && arr[0] == 9.0 );
	}
}

BOOST_AUTO_TEST_CASE(initializer_list_2d) {
#if !defined(__circle_build__) || (__circle_build__ > 200 )  // crashes circle 187-200 in docker
	{
		multi::static_array const arr({
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
		});
		BOOST_TEST_REQUIRE( multi::rank<decltype(arr)>{} == 2 );
		BOOST_TEST_REQUIRE( num_elements(arr) == 6 );
	}
	{
		multi::array const arr({
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
		});
		BOOST_TEST_REQUIRE( multi::rank<decltype(arr)>{} == 2 );
		BOOST_TEST_REQUIRE( num_elements(arr) == 6 );
	}
#endif
}
#endif

BOOST_AUTO_TEST_CASE(partially_formed) {
	multi::array<double, 2> arr1({10, 10}, double{});
	multi::array<double, 2> arr2({10, 10}, {});
	multi::array<double, 2> arr3({10, 10}, 0.0);

	BOOST_REQUIRE( arr1[0][0] == 0.0);
	BOOST_REQUIRE( arr2[0][0] == 0.0);
	BOOST_REQUIRE( arr3[0][0] == 0.0);
}

BOOST_AUTO_TEST_CASE(partially_formed_int_1) {
	multi::array<int, 2> arr1({10, 10}, static_cast<int>(1U));
	multi::array<int, 2> arr2({10, 10}, {1});
	multi::array<int, 2> arr3({10, 10}, 1);

	BOOST_REQUIRE( arr1[0][0] == 1);
	BOOST_REQUIRE( arr2[0][0] == 1);
	BOOST_REQUIRE( arr3[0][0] == 1);
}

BOOST_AUTO_TEST_CASE(partially_formed_int_0) {
	multi::array<int, 2> arr1({10, 10}, int{});
	multi::array<int, 2> arr2({10, 10}, {});
	multi::array<int, 2> arr3({10, 10}, 0);

	BOOST_REQUIRE( arr1[0][0] == 0);
	BOOST_REQUIRE( arr2[0][0] == 0);
	BOOST_REQUIRE( arr3[0][0] == 0);
}

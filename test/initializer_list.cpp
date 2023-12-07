// Copyright 2019-2023 Alfredo A. Correa

#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include <array>
#include <complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d) {
	{
		std::vector<double> const vec = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_REQUIRE( vec[1] == 2. );
	}
	{
		multi::static_array<double, 1> arr = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 5.6 );
	}
	{
		multi::static_array<double, 1> const arr = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 5.6 );
	}
	{
		auto const il = {1.2, 3.4, 5.6};

		multi::static_array<double, 1> const arr(il);
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == il.begin()[2] );
	}
	{
		auto const il = {1.2, 3.4, 5.6};

		multi::static_array<double, 1> const arr(begin(il), end(il));
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == il.begin()[2] );
	}
	{
		multi::static_array<double, 1> const arr = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 5.6 );
		BOOST_REQUIRE(( arr == multi::static_array<double, 1>{1.2, 3.4, 5.6} ));
		BOOST_REQUIRE(( arr == decltype(arr){1.2, 3.4, 5.6} ));
	}
	{
		auto const values = {1.2, 3.4, 5.6};

		multi::array<double, 1> const arr(values.begin(), values.end());
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 5.6 );
	}
	{
		multi::array<double, 1> const arr = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 5.6 );
		BOOST_REQUIRE(( arr == multi::array<double, 1>{1.2, 3.4, 5.6} ));
		BOOST_REQUIRE(( arr == decltype(arr){1.2, 3.4, 5.6} ));
		BOOST_REQUIRE(( arr == decltype(arr)::decay_type({1.2, 3.4, 5.6}) ));
	}
	{
		std::array<double, 3> const stdarr = {
			{1.1, 2.2, 3.3}
                                                                                                                                                                                                      };
		using multi::num_elements;
		BOOST_REQUIRE( num_elements(stdarr) == 3 );

		using std::begin;
		using std::end;
		multi::static_array<double, 1> const arr(begin(stdarr), end(stdarr));
		BOOST_REQUIRE( size(arr) == 3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d_ctad) {
#if defined(__cpp_deduction_guides) and not defined(__NVCC__)
	{
#if not defined(__circle_build__)  // crashes circle 198
		multi::static_array const arr = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 5.6 );
		BOOST_REQUIRE(( arr == multi::static_array{1.2, 3.4, 5.6} ));
#endif
	}
	{
		multi::array arr({1.2, 3.4, 5.6});
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[2] == 5.6 );
		BOOST_REQUIRE(( arr == multi::array({1.2, 3.4, 5.6}) ));
	}
#endif
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_array) {
#if defined(__INTEL_COMPILER) or (defined(__clang__) and (__clang_major__ >= 10))  // doesn't work on gcc
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
		multi::static_array<double, 1> const arr = {1.1, 2.2, 3.3};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( arr[1] == 2.2 );
	}
	{
#if defined(__cpp_deduction_guides) and not defined(__NVCC__)
//      multi::array arr = {{1.1, 2.2, 3.3}};
//      static_assert( decltype(arr)::dimensionality == 1 , "!");
//      BOOST_REQUIRE( size(arr)==3 and arr[1] == 2.2 );
#endif
	}
	{
		std::array<double, 3> stdarr = {{1.1, 2.2, 3.3}};
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
		multi::array<double, 2> const arr = {
			{ 1.2,  2.4, 3.6},
			{11.2, 34.4, 5.6},
			{15.2, 32.4, 5.6},
		};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( size(arr) == 3 and size(arr[0]) == 3 );
		BOOST_REQUIRE( arr[1][1] == 34.4 );
	}
	{
		multi::array<double, 2> arr = {
			{ 1.2,  2.4, 3.6},
			{11.2, 34.4, 5.6},
			{15.2, 32.4, 5.6},
		};
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( size(arr) == 3 and size(arr[0]) == 3 );
		BOOST_REQUIRE( arr[1][1] == 34.4 );
		arr = {
			{ 0.0,  1.0, 02.0},
			{10.0, 11.0, 12.0},
			{20.0, 21.0, 22.0},
		};
		BOOST_REQUIRE( arr[1][2] == 12.0 );
	}
	{
		multi::array<double, 1> vec;
		vec = {4.0, 5.5};
		BOOST_REQUIRE( size(vec) == 2 );
		BOOST_REQUIRE( vec[1] == 5.5 );
	}
	{
		std::array<std::array<double, 2>, 3> const nested = {{{{1.2, 2.4}}, {{11.2, 34.4}}, {{15.2, 32.4}}}};

		using std::begin;
		using std::end;
		multi::static_array<double, 2> arr(begin(nested), end(nested));

		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( size(arr[0]) == 2 );
		BOOST_REQUIRE( arr[1][0] == 11.2 );
	}
	{
		std::array<std::array<double, 2>, 3> const nested = {{{{1.2, 2.4}}, {{11.2, 34.4}}, {{15.2, 32.4}}}};
		multi::static_array<double, 2> const arr(std::begin(nested), std::end(nested));

		BOOST_REQUIRE((
			arr == multi::array<double, 2> {{
				{{ 1.2,  2.4}},
				{{11.2, 34.4}},
				{{15.2, 32.4}}
			}}
		));
		BOOST_REQUIRE(not( arr != multi::array<double, 2>{
				{ 1.2,  2.4},
				{11.2, 34.4},
				{15.2, 32.4}
			}
		));
		BOOST_REQUIRE((
			arr == decltype(arr){
				{ 1.2,  2.4},
				{11.2, 34.4},
				{15.2, 32.4}
			}
		));
	}
	{
		std::array<std::array<double, 2>, 3> nested = {
			{{{1.0, 2.0}},
			 {{2.0, 4.0}},
			 {{3.0, 6.0}}},
		};
		multi::array<double, 2> arr(begin(nested), end(nested));
		BOOST_REQUIRE( num_elements(arr) == 6 and arr[2][1] == 6.0 );
	}
	{
		using complex = std::complex<double>;
		complex const            I{0., 1.};  // NOLINT(readability-identifier-length) imaginary unit
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
	multi::array<double, 3> const arr = {
		{  {1.2, 0.0},  {2.4, 1.0}},
		{ {11.2, 3.0}, {34.4, 4.0}},
		{{15.2, 99.0}, {32.4, 2.0}},
	};
	BOOST_REQUIRE( arr[1][1][0] == 34.4 and arr[1][1][1] == 4.0   );
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d_string) {
	{
		using std::string;
		multi::array<string, 3> B3 = {
			{{"000", "001", "002"}, {"010", "011", "012"}},  // NOLINT(fuchsia-default-arguments-calls)
			{{"100", "101", "102"}, {"110", "111", "112"}},  // NOLINT(fuchsia-default-arguments-calls)
		};
		BOOST_REQUIRE( num_elements(B3)==12 and B3[1][0][1] == "101" );
	}
}

#if defined(__cpp_deduction_guides) and not defined(__NVCC__)
BOOST_AUTO_TEST_CASE(initializer_list_1d_static) {
	{
	#if not defined(__circle_build__)
		multi::static_array arr({1.0, 2.0, 3.0});
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( size(arr) == 3 and num_elements(arr) == 3 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 and num_elements(arr)==3 and arr[1] == 2.0 );
		static_assert(typename decltype(arr)::rank{} == 1);
	#endif
	}
}

BOOST_AUTO_TEST_CASE(initializer_list_1d) {
	{
		multi::array arr({1.0, 2.0, 3.0});
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( size(arr) == 3 and num_elements(arr) == 3 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 and num_elements(arr)==3 and arr[1] == 2.0 );
		static_assert(typename decltype(arr)::rank{} == 1);
	}
	{
		multi::array arr({1.0, 2.0});
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( size(arr) == 2 and num_elements(arr) == 2 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 and num_elements(arr) == 2 and arr[1] == 2.0 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 );
	}
	{
		multi::array arr({0, 2});  //  multi::array arr = {0, 2}; not working with CTAD
		static_assert(std::is_same_v<decltype(arr)::element_type, int>);
		BOOST_REQUIRE( size(arr) == 2 and num_elements(arr) == 2 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 and num_elements(arr) == 2 and arr[1] == 2.0 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 );
	}
	{
		multi::array arr({9.0});  // multi::array arr = {9.0}; not working with CTAD
		static_assert(std::is_same_v<decltype(arr)::element_type, double>);
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 and num_elements(arr)==1 and arr[0]==9. );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{}==1 );
	}
	{
		multi::array arr({9});  // multi::array arr = {9}; not working with CTAD
		static_assert(std::is_same_v<decltype(arr)::element_type, int>);
		BOOST_REQUIRE( size(arr) == 1 and num_elements(arr) == 1 );
		BOOST_REQUIRE( multi::rank<decltype(arr)>{} == 1 );
		BOOST_REQUIRE( num_elements(arr) == 1 and arr[0] == 9.0 );
	}
}

BOOST_AUTO_TEST_CASE(initializer_list_2d) {
	{
	#if not defined(__circle_build__)
		multi::static_array const arr({
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
		});
		BOOST_TEST_REQUIRE( multi::rank<decltype(arr)>{} == 2 );
		BOOST_TEST_REQUIRE( num_elements(arr) == 6 );
	#endif
	}
	{
	#if not defined(__circle_build__)
		multi::array const arr({
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
		});
		BOOST_TEST_REQUIRE( multi::rank<decltype(arr)>{} == 2 );
		BOOST_TEST_REQUIRE( num_elements(arr) == 6 );
	#endif
	}
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
	multi::array<int, 2> arr1({10, 10}, int{1});
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

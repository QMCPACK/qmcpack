// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, static_array, num_elements

// IWYU pragma: no_include <algorithm>  // for copy  // bug in iwyu 14.0.6? with GNU stdlib
#include <array>             // for array
#include <complex>           // for operator*, operator+, complex
#include <initializer_list>  // for initializer_list, begin, end
#include <iterator>          // for size, begin, end
#include <string>            // for basic_string, allocator, char_tr...
#include <type_traits>       // for is_same_v
#include <vector>            // for vector

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d) {
		{
			std::vector<int> const vec = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)
			BOOST_TEST( vec[1] == 20 );
		}
		{
			multi::static_array<int, 1> arr = {12, 34, 56};
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
		}
		{
			multi::static_array<int, 1> const arr = {12, 34, 56};
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
		}
		{
			auto const il = {12, 34, 56};

			multi::static_array<int, 1> const arr(il);
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		}
		{
			auto const il = {12, 34, 56};

			multi::static_array<int, 1> const arr(begin(il), end(il));
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		}
		{
			multi::static_array<int, 1> const arr = {12, 34, 56};
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
			BOOST_TEST(( arr == multi::static_array<int, 1>{12, 34, 56} ));
			BOOST_TEST(( arr == decltype(arr){12, 34, 56} ));
		}
		{
			auto const values = {12, 34, 56};

			multi::array<int, 1> const arr(values.begin(), values.end());
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
		}
		{
			multi::array<int, 1> const arr = {12, 34, 56};

			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );

			BOOST_TEST(( arr == multi::array<int, 1>{12, 34, 56} ));
			BOOST_TEST(( arr == decltype(arr){12, 34, 56} ));
			BOOST_TEST(( arr == decltype(arr)::decay_type({12, 34, 56}) ));
		}
		{
			std::array<int, 3> const stdarr = {
				{11, 22, 33},
			};
			using multi::num_elements;
			BOOST_TEST( num_elements(stdarr) == 3 );

			using std::begin;
			using std::end;
			multi::static_array<double, 1> const arr(begin(stdarr), end(stdarr));
			BOOST_TEST( size(arr) == 3 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d_ctad) {
#if defined(__cpp_deduction_guides) && !defined(__NVCC__)
		// #if __cplusplus >= 202002L
		// static constexpr auto f = []
		// {
		//  multi::array<int, 1> arr(3);
		//  arr[0] = 12; arr[1] = 34; arr[2] = 56;  // TODO(correaa) getting "assignment to object outside its lifetime is not allowed in a constant expression"
		//  return (arr.size() == 3);
		// }();
		// static_assert(f);
		// #endif

		{
			multi::static_array const arr = {12, 34, 56};
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
			BOOST_TEST(( arr == multi::static_array{12, 34, 56} ));
		}
		{
			multi::array<int, 1> arr(std::initializer_list<int>{12, 34, 56});
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
			BOOST_TEST(( arr == multi::array<int, 1>(std::initializer_list<int>{12, 34, 56}) ));
		}
		{
			multi::array<int, 1> arr({12, 34, 56});
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
			BOOST_TEST(( arr == multi::array<int, 1>({12, 34, 56}) ));
		}
		#if !defined(__GNUC__) || (__GNUC__ < 14)  // workaround bug in gcc 14.2
		{
			multi::array arr({12, 34, 56});
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
			BOOST_TEST(( arr == multi::array({12, 34, 56}) ));
		}
		#endif
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
			BOOST_TEST( arr.size() == 11 );
			BOOST_TEST( arr[9] == 9.0 );
		}
#endif
	}

	BOOST_AUTO_TEST_CASE(multi_initialize_from_carray_1d) {
		{
			multi::static_array<int, 1> const arr = {11, 22, 33};
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[1] == 22 );
		}
		{
#if defined(__cpp_deduction_guides) && !defined(__NVCC__)
//      multi::array arr = {{1.1, 2.2, 3.3}};
//      static_assert( decltype(arr)::dimensionality == 1 , "!");
//      BOOST_TEST( size(arr)==3 && arr[1] == 2.2 );
#endif
		}
		{
			std::array<double, 3> stdarr = {
				{1.1, 2.2, 3.3}
			};
			multi::array<double, 1> const arr(begin(stdarr), end(stdarr));
			BOOST_TEST(( arr == decltype(arr){1.1, 2.2, 3.3} ));
		}
	}

	BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_2d) {
		{
			multi::static_array<double, 2> const arr = {
				{ 1.2,  2.4, 3.6, 8.9},
				{11.2, 34.4, 5.6, 1.1},
				{15.2, 32.4, 5.6, 3.4},
			};
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( size(arr[0]) == 4 );
			BOOST_TEST(( arr == decltype(arr){
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

			multi::static_array<int, 2> arr(begin(nested), end(nested));

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
			std::array<std::array<int, 2>, 3> const nested = {
				{
                 {{10, 20}},
                 {{20, 40}},
                 {{30, 60}},
				 }
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
			BOOST_TEST( arr[1][1] == 4.0 + 2.0*I );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_tests_static_array_initializer_list) {
		multi::static_array<std::complex<double>, 2> SA = {
			{{1.0, 0.0}, {2.0, 0.0}},
			{{3.0, 0.0}, {4.0, 0.0}},
		};
		BOOST_TEST( SA[1][1] == 4.0 );
	}

	BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d) {
		multi::array<int, 3> const arr = {
			{ {12, 100},  {24, 10}},
			{ {112, 30}, {344, 40}},
			{{152, 990}, {324, 20}},
		};
		BOOST_TEST( arr[1][1][0] == 344 );
		BOOST_TEST( arr[1][1][1] ==  40 );
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

			BOOST_TEST( num_elements(B3) == 12 );
			BOOST_TEST( B3[1][0][1] == "101" );
		}
	}

#if defined(__cpp_deduction_guides) && !defined(__NVCC__)
	BOOST_AUTO_TEST_CASE(initializer_list_1d_static) {
		multi::static_array arr({10, 20, 30});

		static_assert(std::is_same_v<decltype(arr)::element_type, int>);

		BOOST_TEST( size(arr) == 3 && num_elements(arr) == 3 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1);
		BOOST_TEST( num_elements(arr) == 3 );
		BOOST_TEST( arr[1] == 20 );

		static_assert(typename decltype(arr)::rank{} == 1);
	}

	#if !defined(__GNUC__) || (__GNUC__ < 14)  // workaround bug in gcc 14.2
	BOOST_AUTO_TEST_CASE(initializer_list_1d_a) {
		multi::array arr({10, 20, 30});

		static_assert(std::is_same_v<decltype(arr)::element_type, int>);

		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( num_elements(arr) == 3 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( num_elements(arr) == 3 );
		BOOST_TEST( arr[1] == 20 );

		static_assert(typename decltype(arr)::rank{} == 1);
	}

	BOOST_AUTO_TEST_CASE(initializer_list_1d_b) {
		multi::array arr({10, 20});
		static_assert(std::is_same_v<decltype(arr)::element_type, int>);

		BOOST_TEST( size(arr) == 2 );
		BOOST_TEST( num_elements(arr) == 2 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( num_elements(arr) == 2 );
		BOOST_TEST( arr[1] == 20 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
	}

	BOOST_AUTO_TEST_CASE(initializer_list_1d_c) {
		multi::array arr({0, 2});  //  multi::array arr = {0, 2}; not working with CTAD

		static_assert(std::is_same_v<decltype(arr)::element_type, int>);

		BOOST_TEST( size(arr) == 2 );
		BOOST_TEST( num_elements(arr) == 2 );
		BOOST_TEST( multi::rank<decltype(arr)>{} == 1 );
		BOOST_TEST( num_elements(arr) == 2 );
		BOOST_TEST( arr[1] == 2 );
		BOOST_TEST( multi::rank<decltype(arr)>{} == 1 );
	}

	BOOST_AUTO_TEST_CASE(initializer_list_1d_d) {
		multi::array arr({90});  // multi::array arr = {90}; not working with CTAD

		static_assert(std::is_same_v<decltype(arr)::element_type, int>);

		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( num_elements(arr) == 1 );
		BOOST_TEST( arr[0] == 90 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
	}

	BOOST_AUTO_TEST_CASE(initializer_list_1d_e) {
		multi::array arr({90});  // multi::array arr = {90}; not working with CTAD

		static_assert(std::is_same_v<decltype(arr)::element_type, int>);

		BOOST_TEST( size(arr) == 1 );
		BOOST_TEST( num_elements(arr) == 1 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( num_elements(arr) == 1 );
		BOOST_TEST( arr[0] == 90 );
	}

	BOOST_AUTO_TEST_CASE(initializer_list_2d) {
		{
			multi::static_array const arr({
				{1.0, 2.0, 3.0},
				{4.0, 5.0, 6.0},
			});
			BOOST_TEST( multi::rank<decltype(arr)>{} == 2 );
			BOOST_TEST( num_elements(arr) == 6 );
		}
		{
			multi::array const arr({
				{1.0, 2.0, 3.0},
				{4.0, 5.0, 6.0},
			});
			BOOST_TEST( multi::rank<decltype(arr)>::value == 2 );
			BOOST_TEST( num_elements(arr) == 6 );
		}
	}
	#endif
#endif

	BOOST_AUTO_TEST_CASE(partially_formed) {
		multi::array<int, 2> arr1({10, 10}, int{});
		multi::array<int, 2> arr2({10, 10}, {});
		multi::array<int, 2> arr3({10, 10}, 0);

		BOOST_TEST( arr1[0][0] == 0);
		BOOST_TEST( arr2[0][0] == 0);
		BOOST_TEST( arr3[0][0] == 0);
	}

	BOOST_AUTO_TEST_CASE(partially_formed_int_1) {
		multi::array<int, 2> arr1({10, 10}, static_cast<int>(1U));
		multi::array<int, 2> arr2({10, 10}, {1});
		multi::array<int, 2> arr3({10, 10}, 1);

		BOOST_TEST( arr1[0][0] == 1);
		BOOST_TEST( arr2[0][0] == 1);
		BOOST_TEST( arr3[0][0] == 1);
	}

	BOOST_AUTO_TEST_CASE(partially_formed_int_0) {
		multi::array<int, 2> arr1({10, 10}, int{});
		multi::array<int, 2> arr2({10, 10}, {});
		multi::array<int, 2> arr3({10, 10}, 0);

		BOOST_TEST( arr1[0][0] == 0);
		BOOST_TEST( arr2[0][0] == 0);
		BOOST_TEST( arr3[0][0] == 0);
	}

	return boost::report_errors();
}

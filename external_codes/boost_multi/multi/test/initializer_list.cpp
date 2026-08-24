// Copyright 2019-2026 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>        // for array, dynamic_array, num_elements
#include <boost/multi/restriction.hpp>  // for array, dynamic_array, num_elements

#include <boost/core/lightweight_test.hpp>

#include <algorithm>         // IWYU pragma: keep  // for copy
#include <array>             // for array
#include <cmath>             // IWYU pragma: keep  // for abs
#include <complex>           // for operator*, operator+, complex
#include <initializer_list>  // for initializer_list, begin, end
#include <iterator>          // for size, begin, end
#include <string>            // for basic_string, allocator, char_tr...
#include <tuple>             // IWYU pragma: keep
#include <type_traits>       // for is_same_v
// IWYU pragma: no_include <utility>           // for declval, forward, move
#include <vector>  // for vector

namespace multi = boost::multi;

namespace boost::multi {

template<class T>
auto operator+(std::initializer_list<T> il) {  // NOLINT(misc-use-anonymous-namespace,misc-use-internal-linkage)
	multi::array<T, 1> ret({static_cast<multi::ssize_t>(il.size())}, T{});
	std::copy(il.begin(), il.end(), ret.begin());
	return ret;
}

template<class T>
auto operator+(std::initializer_list<std::initializer_list<T>> il) {  // NOLINT(misc-use-anonymous-namespace,misc-use-internal-linkage)
	auto const size2 = il.size() == 0 ? 0 : std::max_element(il.begin(), il.end(), [](auto const& a, auto const& b) { return a.size() < b.size(); })->size();

	multi::array<T, 2> ret({static_cast<multi::ssize_t>(il.size()), static_cast<multi::ssize_t>(size2)}, T{});
	std::copy(il.begin(), il.end(), ret.begin());
	return ret;
}

}  // end namespace boost::multi

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d)
	{
		std::vector<int> const vec = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( vec[1] == 20 );
	}
	{
		multi::dynamic_array<int, 1> arr = {12, 34, 56};
		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( arr[2] == 56 );
	}
	{
		multi::dynamic_array<int, 1> const arr = {12, 34, 56};
		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( arr[2] == 56 );
	}
	{
		auto const il = {12, 34, 56};

		multi::dynamic_array<int, 1> const arr(il);
		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( arr[2] == 56 );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
	{
		auto const il = {12, 34, 56};

		multi::dynamic_array<int, 1> const arr(begin(il), end(il));
		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( arr[2] == 56 );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
	{
		multi::dynamic_array<int, 1> const arr = {12, 34, 56};
		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( arr[2] == 56 );
		BOOST_TEST(( arr == multi::dynamic_array<int, 1>{12, 34, 56} ));
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
	}
	{
		std::array<int, 3> const stdarr = {
			{11, 22, 33},
		};
		using multi::num_elements;
		BOOST_TEST( num_elements(stdarr) == 3 );

		using std::begin;
		using std::end;
		multi::dynamic_array<double, 1> const arr(begin(stdarr), end(stdarr));
		BOOST_TEST( size(arr) == 3 );
	}

	// BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d_ctad)
	{
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
			// multi::dynamic_array const arr = {12, 34, 56};
			multi::dynamic_array<int, 1> const arr = {12, 34, 56};

			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
			BOOST_TEST(( arr == multi::dynamic_array<int, 1>{12, 34, 56} ));
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
			// multi::array arr({12, 34, 56});
			multi::array<int, 1> arr({12, 34, 56});
			BOOST_TEST( size(arr) == 3 );
			BOOST_TEST( arr[2] == 56 );
			BOOST_TEST(( arr == multi::array<int, 1>({12, 34, 56}) ));
		}
#endif
#endif
	}

	// BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_array)
	{
#if (defined(__INTEL_COMPILER) || (defined(__clang__) && (__clang_major__ >= 10))) && !defined(__circle_build__)  // doesn't work on gcc; circle rejects designator + trailing positional initializers ("list comprehension disabled")
		{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-designator"
			//      double const a[] = { [8] = 8.0, 9.0, 10.0 };
			std::array<double, 11> const stdarr = {
				{[8] = 8.0, 9.0, 10.0},
			};
#pragma GCC diagnostic pop
			multi::array<double, 1> arr(stdarr);
			BOOST_TEST( arr.size() == 11 );
			BOOST_TEST( arr[9] == 9.0 );
		}
#endif
	}

	// BOOST_AUTO_TEST_CASE(multi_initialize_from_carray_1d)
	{
		multi::dynamic_array<int, 1> const arr = {11, 22, 33};
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

	// BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_2d)
	{
		multi::dynamic_array<double, 2> const arr = {
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
		BOOST_TEST( arr.size() == 3 );
		BOOST_TEST( arr[0].size() == 3 );
		BOOST_TEST( arr[1][1] == 344 );

		using multi::operator+;  // cppcheck-suppress [constStatement];

		auto arr2 = operator+({
			{ 12,  24, 36},
			{112, 344, 56},
			{152, 324, 56},
		});

		BOOST_TEST( arr2 == arr );
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

		multi::dynamic_array<int, 2> arr(begin(nested), end(nested));

		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( size(arr[0]) == 2 );
		BOOST_TEST( arr[1][0] == 112 );
	}
	{
		std::array<std::array<int, 2>, 3> const nested = {
			{{{12, 24}}, {{112, 344}}, {{152, 324}}}
		};
		multi::dynamic_array<int, 2> const arr(std::begin(nested), std::end(nested));

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
		BOOST_TEST( arr.num_elements() == 6 );
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

	// BOOST_AUTO_TEST_CASE(multi_tests_dynamic_array_initializer_list)
	{
		multi::dynamic_array<std::complex<double>, 2> SA = {
			{{1.0, 0.0}, {2.0, 0.0}},
			{{3.0, 0.0}, {4.0, 0.0}},
		};
		BOOST_TEST( SA[1][1] == 4.0 );
	}

	// BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d)
	{
		multi::array<int, 3> const arr = {
			{ {12, 100},  {24, 10}},
			{ {112, 30}, {344, 40}},
			{{152, 990}, {324, 20}},
		};
		BOOST_TEST( arr[1][1][0] == 344 );
		BOOST_TEST( arr[1][1][1] ==  40 );
	}

	// BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d_string)
	{
		{
			using std::string;

			// NOLINTBEGIN(fuchsia-default-arguments-calls)
			multi::array<string, 3> B3 = {
				{{"000", "001", "002"}, {"010", "011", "012"}},
				{{"100", "101", "102"}, {"110", "111", "112"}},
			};
			// NOLINTEND(fuchsia-default-arguments-calls)

			BOOST_TEST( B3.num_elements() == 12 );
			BOOST_TEST( B3[1][0][1] == "101" );
		}
	}

#if defined(__cpp_deduction_guides) && !defined(__NVCC__)
	// BOOST_AUTO_TEST_CASE(initializer_list_1d_static)
	{
		// multi::dynamic_array arr({10, 20, 30});
		multi::dynamic_array<int, 1> arr({10, 20, 30});

		static_assert(std::is_same_v<decltype(arr)::element, int>);

		BOOST_TEST( size(arr) == 3 && arr.num_elements() == 3 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1);
		BOOST_TEST( arr.num_elements() == 3 );
		BOOST_TEST( arr[1] == 20 );

		static_assert(decltype(arr)::rank{} == 1);
	}

#if !defined(__GNUC__) || (__GNUC__ < 14)  // workaround bug in gcc 14.2
	// BOOST_AUTO_TEST_CASE(initializer_list_1d_a)
	{
		multi::array<int, 1> arr({10, 20, 30});

		static_assert(std::is_same_v<decltype(arr)::element, int>);

		BOOST_TEST( size(arr) == 3 );
		BOOST_TEST( arr.num_elements() == 3 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( arr.num_elements() == 3 );
		BOOST_TEST( arr[1] == 20 );

		static_assert(decltype(arr)::rank{} == 1);
	}

	// BOOST_AUTO_TEST_CASE(initializer_list_1d_b)
	{
		multi::array<int, 1> arr({10, 20});
		static_assert(std::is_same_v<decltype(arr)::element, int>);

		BOOST_TEST( arr.size() == 2 );
		BOOST_TEST( arr.num_elements() == 2 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( arr.num_elements() == 2 );
		BOOST_TEST( arr[1] == 20 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
	}

	// BOOST_AUTO_TEST_CASE(initializer_list_1d_c)
	{
		// multi::array arr({0, 2});  //  multi::array arr = {0, 2}; not working with CTAD
		multi::array<int, 1> arr({0, 2});

		static_assert(std::is_same_v<decltype(arr)::element, int>);

		BOOST_TEST( size(arr) == 2 );
		BOOST_TEST( arr.num_elements() == 2 );
		BOOST_TEST( multi::rank<decltype(arr)>{} == 1 );
		BOOST_TEST( arr.num_elements() == 2 );
		BOOST_TEST( arr[1] == 2 );
		BOOST_TEST( multi::rank<decltype(arr)>{} == 1 );
	}

	// BOOST_AUTO_TEST_CASE(initializer_list_1d_d)
	{
		// multi::array arr({90});  // multi::array arr = {90}; not working with CTAD
		multi::array<int, 1> arr({90});  // multi::array arr = {90}; not working with CTAD

		static_assert(std::is_same_v<decltype(arr)::element, int>);

		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( arr.num_elements() == 1 );
		BOOST_TEST( arr[0] == 90 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
	}

	// BOOST_AUTO_TEST_CASE(initializer_list_1d_e)
	{
		// multi::array arr({90});  // multi::array arr = {90}; not working with CTAD
		multi::array<int, 1> arr({90});  // multi::array arr = {90}; not working with CTAD

		static_assert(std::is_same_v<decltype(arr)::element, int>);

		BOOST_TEST( size(arr) == 1 );
		BOOST_TEST( arr.num_elements() == 1 );
		BOOST_TEST( multi::rank<decltype(arr)>::value == 1 );
		BOOST_TEST( arr.num_elements() == 1 );
		BOOST_TEST( arr[0] == 90 );
	}

	// BOOST_AUTO_TEST_CASE(initializer_list_2d)
	{
		{
			multi::dynamic_array<double, 2> const arr = {
				{1.0, 2.0, 3.0},
				{4.0, 5.0, 6.0},
			};

			BOOST_TEST( decltype(arr)::dimensionality == 2 );
			// BOOST_TEST( multi::rank<decltype(arr)>{} == 2 );
			BOOST_TEST( arr.num_elements() == 6 );
		}
	}
#endif
#endif

	// BOOST_AUTO_TEST_CASE(partially_formed)
	{
		multi::array<int, 2> arr1({10, 10}, int{});
		multi::array<int, 2> arr2({10, 10}, {});
		multi::array<int, 2> arr3({10, 10}, 0);

		BOOST_TEST( arr1[0][0] == 0);
		BOOST_TEST( arr2[0][0] == 0);
		BOOST_TEST( arr3[0][0] == 0);
	}

	// BOOST_AUTO_TEST_CASE(partially_formed_int_1)
	{
		multi::array<int, 2> arr1({10, 10}, static_cast<int>(1U));
		multi::array<int, 2> arr2({10, 10}, {1});
		multi::array<int, 2> arr3({10, 10}, 1);

		BOOST_TEST( arr1[0][0] == 1);
		BOOST_TEST( arr2[0][0] == 1);
		BOOST_TEST( arr3[0][0] == 1);
	}

	// BOOST_AUTO_TEST_CASE(partially_formed_int_0)
	{
		multi::array<int, 2> arr1({10, 10}, int{});
		multi::array<int, 2> arr2({10, 10}, {});
		multi::array<int, 2> arr3({10, 10}, 0);

		BOOST_TEST( arr1[0][0] == 0);
		BOOST_TEST( arr2[0][0] == 0);
		BOOST_TEST( arr3[0][0] == 0);
	}

#ifdef __cpp_deduction_guides
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif
	{
		std::initializer_list<int> const il = {1, 2, 3};

		BOOST_TEST(*multi::base(il) == 1);

		multi::const_subarray<int, 1> const csarr(il);

		BOOST_TEST( csarr.size() == 3 );
		BOOST_TEST( csarr.num_elements() == 3 );

		BOOST_TEST( csarr[1] == il.begin()[1] );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

		multi::const_subarray const csarr2(il);

		BOOST_TEST( csarr2 == csarr );
	}
#ifdef __clang__
#pragma clang diagnostic pop
#endif
	{
		BOOST_TEST( multi::const_subarray({1, 2, 3})[1] == 2 );
	}
#endif
	{
		multi::layout_t<2> const lyt(multi::extents_t<2>(3, 2));
		BOOST_TEST( lyt.num_elements() == 6 );
	}
	// {
	// 	std::initializer_list<std::initializer_list<int>> const il = {
	// 		{11}, {22}, {33},
	// 	};

	// 	auto bas = multi::base(il);

	// 	BOOST_TEST(*bas == 11);
	// 	BOOST_TEST( multi::layout(il).size() == 3 );

	// 	multi::const_subarray<int, 2> const csarr(il);

	// 	BOOST_TEST( csarr.size() == 3 );
	// 	BOOST_TEST( csarr[0].size() == 1 );
	// 	BOOST_TEST( csarr[1].size() == 1 );
	// 	BOOST_TEST( csarr[2].size() == 1 );

	// 	BOOST_TEST( csarr[0][0] == 11 );
	// 	BOOST_TEST( csarr[1][0] == 22 );
	// 	std::cout << "stride = " << csarr.stride() << std::endl;
	// 	std::cout << "list = " << " *" << bas[0] << ' ' << bas[csarr.stride()] << ' ' << bas[2*csarr.stride()] << " " << bas[2*csarr.stride()] << std::endl;
	// 	if(csarr[2][0] != 33) {
	// 		for(int i = -32; i != 64; ++i) {
	// 			std::cout << "bas[" << i << "] = " << bas[i] << std::endl;
	// 		}
	// 	}
	// 	BOOST_TEST( csarr[2][0] == 33 );
	// }
	{
		// std::initializer_list<std::initializer_list<int>> const il = {
		// 	{1, 2, 3},
		// 	{4, 5, 6}
		// };

		// static_assert(std::is_same_v<multi::element_t<std::decay_t<decltype(il)>>, int>);

		// BOOST_TEST(*multi::base(il) == 1);

		// auto const il_lyt = multi::layout(il);

		// auto [s1, s2] = il_lyt.strides();

		// BOOST_TEST( std::abs(s1) >= 3 );
		// BOOST_TEST( s2 == 1 );

		// BOOST_TEST( il_lyt.size() == 2);
		// BOOST_TEST( il_lyt.num_elements() == 6 );

		// multi::const_subarray<int, 2> const csarr(il);

		// using std::get;
		// BOOST_TEST( get<0>(csarr.sizes()) == 2 );
		// BOOST_TEST( get<1>(csarr.sizes()) == 3 );

		// BOOST_TEST( csarr[0][0] == 1);
		// BOOST_TEST( csarr[0][1] == 2);
		// BOOST_TEST( csarr[0][2] == 3);

		// BOOST_TEST( csarr[1][0] == 4);
		// BOOST_TEST( csarr[1][1] == 5);
		// BOOST_TEST( csarr[1][2] == 6);

		// multi::array<int, 2> const arr = csarr;

		// BOOST_TEST( arr == csarr );

		// multi::array<int, 2> const arr2{multi::const_subarray<int, 2>(il)};

		// BOOST_TEST( arr == arr2 );

		// multi::array<int, 2> const arr3 = multi::const_subarray<int, 2>(il);

		// BOOST_TEST( arr == arr3 );

		// multi::const_subarray const csarr2(il);

		// BOOST_TEST( csarr2 == csarr );

		// multi::array<int, 2> const arr4 = multi::const_subarray(il);

		// BOOST_TEST( arr == arr4 );

		// multi::array const arr5 = multi::const_subarray(il);

		// BOOST_TEST( arr == arr5 );

		// auto const arr6 = +multi::const_subarray(il);

		// BOOST_TEST( arr == arr6 );

		// using multi::operator+;
		// auto const arr7 = +il;

		// BOOST_TEST( arr == arr7 );

		// // +{...} doesn't compile
		// auto const arr8 = operator+({
		// 	{1, 2, 3},
		// 	{4, 5, 6}
		// });

		// BOOST_TEST( arr == arr8 );
	}
	{
		std::initializer_list<int> const il = {};
		// BOOST_TEST(il.begin() == nullptr);  // depends on the implementation
		BOOST_TEST(il.size() == 0);  // depends on the implementation
	}
	{
		std::initializer_list<std::initializer_list<int>> const il = {
			{1, 2, 3}
		};

		BOOST_TEST( il.size() == 1 );
		BOOST_TEST( il.begin()->size() == 3 );

		BOOST_TEST(multi::base(il) != nullptr);
		auto const* bas = multi::base(il);
		if(bas != nullptr) {
			BOOST_TEST(*bas == 1);
		}
		// std::cout << "size = " << multi::layout(il).size() << std::endl;
		// BOOST_TEST( multi::layout(il).size() == 1 );

		// multi::const_subarray<int, 2> const csarr(il);

		// BOOST_TEST( csarr.size() == 1 );
		// BOOST_TEST( csarr.num_elements() == 3 );

		// BOOST_TEST( csarr[0].size() == 3 );
		// BOOST_TEST( csarr[0][0] == 1 );
		// BOOST_TEST( csarr[0][1] == 2 );
		// BOOST_TEST( csarr[0][2] == 3 );
	}
	{
		std::initializer_list<std::initializer_list<int>> const il = {
			{1, 2, 3},
			{4, 5, 6}
		};

		BOOST_TEST( il.size() == 2 );
		BOOST_TEST( il.begin()->size() == 3 );

		multi::array<int, 2> arr = multi::detail::make_restriction(il);

		// multi::const_subarray<int, 2> const csarr(il);

		BOOST_TEST( arr.size() == 2 );
		BOOST_TEST( arr.num_elements() == 6 );

		BOOST_TEST( arr[0].size() == 3 );

		BOOST_TEST( arr[0][0] == 1 );
		BOOST_TEST( arr[0][1] == 2 );
		BOOST_TEST( arr[0][2] == 3 );

		BOOST_TEST( arr[1][0] == 4 );
		BOOST_TEST( arr[1][1] == 5 );
		BOOST_TEST( arr[1][2] == 6 );
	}
	{
		std::initializer_list<std::initializer_list<int>> const il = {
			{11},
			{22},
			{33},
		};

		multi::array<int, 2> const arr = multi::detail::make_restriction(il);

		BOOST_TEST( arr.size() == 3 );
		BOOST_TEST( arr[0].size() == 1 );
		BOOST_TEST( arr[1].size() == 1 );
		BOOST_TEST( arr[2].size() == 1 );

		BOOST_TEST( arr[0][0] == 11 );
		BOOST_TEST( arr[1][0] == 22 );

		BOOST_TEST( arr[0][0] == 11 );
		BOOST_TEST( arr[1][0] == 22 );
		BOOST_TEST( arr[2][0] == 33 );
	}

	// NOLINTBEGIN(readability-identifier-length,misc-const-correctness)
	// Andrzej examples
	{
		{
			// multi::array A = {2, 3};  // doesn't compile anymore (Jan 19 2026)
			// multi::array B = {{2, 3}, {4, 5}};  // doesn't compile anymore (Jan 19 2026)
		}
		{
			multi::array<int, 2> A({2, 3});
			multi::array<int, 2> B{
				{2, 3}
			};

			using std::get;
			BOOST_TEST( get<0>(A.sizes()) == 2 );
			BOOST_TEST( get<1>(A.sizes()) == 3 );

			BOOST_TEST( get<0>(B.sizes()) == 1 );
			BOOST_TEST( get<1>(B.sizes()) == 2 );
			BOOST_TEST( B[0][0] == 2 );
			BOOST_TEST( B[0][1] == 3 );
		}
		{
			// multi::array A1 ({2, 3});  // doesn't compile as of Jan 29
			// multi::array B2 {{2, 3}};  // doesn't compile as of Jan 29
		}
		{
			multi::array<int, 1> C(2);
			multi::array<int, 1> D{2};

			BOOST_TEST( C.size() == 2);
			BOOST_TEST( D.size() == 1);
			BOOST_TEST( D[0] == 2 );
		}
		{
			multi::array<int, 2> A1({2, 3});  // argument interpreted as extents
			// multi::array         A2 ( {2, 3} );  // doesn't compile as of Jan 29

			using std::get;
			BOOST_TEST( get<0>(A1.sizes()) == 2 );
			BOOST_TEST( get<1>(A1.sizes()) == 3 );
		}
		{
			multi::array<int, 1> A1({3}, 11);
			multi::array<int, 1> A2({3});

			BOOST_TEST( A1.size() == 3 );
			BOOST_TEST( A1[0] == 11 );

			BOOST_TEST( A2.size() == 1 );
			BOOST_TEST( A2[0] == 3 );
		}
	}
	// NOLINTEND(readability-identifier-length,misc-const-correctness)

	// #ifndef __circle_build__
	// 	{
	// 		BOOST_TEST( (multi::initializer_array<int, 2>{
	// 			{1, 2, 3},
	// 			{4, 5, 6}
	// 		})[1][1] == 5 );
	// 	}
	// 	{
	// 		BOOST_TEST( (multi::initializer_array<int, 1>{1, 2, 3})[1] == 2 );
	// 	}
	// #endif
	{
		multi::array<int, 1> numbers(3);
		numbers = {1, 2, 3};
	}

	return boost::report_errors();
}

// Copyright 2021-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>  // for array, layout_t, subarray, sizes

#include <array>    // for array
#include <numeric>  // for iota
#if(__cplusplus >= 202002L)
	#include <ranges>  // IWYU pragma: keep
#endif
#include <type_traits>  // for is_assignable_v

namespace multi = boost::multi;

#if(__cplusplus >= 202002L)
	#if defined(__cpp_lib_ranges_repeat) && (__cpp_lib_ranges_repeat >= 202207L)

template<class X1D, class Y1D>
auto meshgrid(X1D const& x, Y1D const& y) {
	return std::pair{x.broadcasted().rotated(), y.broadcasted()};
}

template<class X1D, class Y1D>
auto meshgrid_copy(X1D const& x, Y1D const& y) {
	auto ret = std::pair{
		multi::array<typename X1D::element_type, 2>({x.size(), y.size()}),
		multi::array<typename Y1D::element_type, 2>(std::views::repeat(y, x.size()))
	};

	std::fill(ret.first.rotated().begin(), ret.first.rotated().end(), x);
	// std::ranges::fill(ret.first.rotated(), x);

	return ret;
}
	#endif
#endif

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
#ifndef _MSC_VER      // msvc 14.40 gets confused with constexpr
	BOOST_AUTO_TEST_CASE(constexpr_carray_rotated_end) {
		constexpr auto test = [] {
			std::array<int,
					   // 9  // produces UB: cannot refer to element 10 of array of 9 elements in a constant expression
					   12  // ok
					   >
				buffer = {
					{0, 1, 2, 3, 4, 5, 6, 7, 8}
            };  // , 10, 11};

			multi::array_ref<int, 2> arr({3, 3}, buffer.data());  // // TODO(correaa) think how to handle references to arrays (UB)

			auto const& arrrot1 = arr.rotated()[1];
			return (arrrot1.end() != arrrot1.begin());
		}();
		BOOST_TEST(test);
	}

	BOOST_AUTO_TEST_CASE(constexpr_carray_rotated_end_interval) {
		constexpr auto test = [] {
			std::array<int,
					   // 9  // produces UB: cannot refer to element 10 of array of 9 elements in a constant expression
					   12  // ok
					   >
				buffer = {
					{0, 1, 2, 3, 4, 5, 6, 7, 8}
            };  // , 10, 11};

			multi::array_ref<int, 2> arr({3, 3}, buffer.data());  // // TODO(correaa) think how to handle references to arrays (UB)

			return arr.rotated()[1]({3, 3}).size() == 0;
		}();
		BOOST_TEST(test);
	}

	BOOST_AUTO_TEST_CASE(constexpr_carray_diagonal_end_2D) {
		constexpr auto test = [] {
			std::array<int,
					   12 + 3  // (3 + 1)*4 - 1
					   >
				buffer = {};

			multi::array_ref<int, 2> arr({3, 4}, buffer.data());  // // TODO(correaa) think how to handle references to arrays (UB)

			return arr.diagonal().end() != arr.diagonal().begin();
		}();
		BOOST_TEST(test);
	}

	BOOST_AUTO_TEST_CASE(constexpr_carray_rotated_end_3D) {
		constexpr auto test = [] {
			std::array<int,
					   60 + 18  // 3*4*5 + (4-1)*5 + (5 - 1) or 60 + 15 + 4 = 79
					   >
				buffer = {};

			multi::array_ref<int, 3> arr({3, 4, 5}, buffer.data());  // // TODO(correaa) think how to handle references to arrays (UB)

			return arr.diagonal().diagonal().end() != arr.diagonal().diagonal().begin();
		}();
		BOOST_TEST(test);
	}

	#if __cplusplus >= 202002L
		#if defined(__GNUC__) && !defined(__clang__) && !defined(__NVCOMPILER) && !defined(__NVCC__)
	BOOST_AUTO_TEST_CASE(constexpr_dynamic_array_rotated_end) {
		constexpr auto test = [] {
			multi::array<int, 2> arr({3, 3}, 99);
			return arr.rotated()[1].end() != arr.rotated()[1].begin();
		}();
		BOOST_TEST(test);
	}
		#endif
	#endif

#endif  // MSVC

	BOOST_AUTO_TEST_CASE(multi_2d_const) {
		multi::array<int, 2> const arr = {
			{10, 20},
			{30, 40},
		};

		BOOST_TEST( arr.rotated()[1][1] == 40 );
		static_assert(!std::is_assignable_v<decltype(arr.rotated()[1][1]), decltype(50)>);
	}

	BOOST_AUTO_TEST_CASE(multi_2d) {
		multi::array<int, 2> arr = {
			{10, 20},
			{30, 40},
		};

		BOOST_TEST( arr.rotated()[1][1] == 40 );

		// what(arr.rotated()[0][0]);

		static_assert(std::is_assignable_v<decltype(arr.rotated()[1][1]), decltype(50)>);

		// arr.rotated()[1][1] = 50;
	}

	BOOST_AUTO_TEST_CASE(multi_rotate_3d) {
		multi::array<double, 3> arr({3, 4, 5});

		using std::get;

		BOOST_TEST( get<0>(arr.sizes()) == 3 );
		BOOST_TEST( get<1>(arr.sizes()) == 4 );
		BOOST_TEST( get<2>(arr.sizes()) == 5 );

		auto&& RA = arr.rotated();
		BOOST_TEST(( sizes(RA) == decltype(RA.sizes()){4, 5, 3} ));
		BOOST_TEST(  &arr[0][1][2] == &RA[1][2][0] );

		auto&& UA = arr.unrotated();
		BOOST_TEST(( sizes(UA) == decltype(sizes(UA)){5, 3, 4} ));
		BOOST_TEST( &arr[0][1][2] == &UA[2][0][1] );

		auto&& RRA = RA.rotated();
		BOOST_TEST(( sizes(RRA) == decltype(sizes(RRA)){5, 3, 4} ));
		BOOST_TEST( &arr[0][1][2] == &RRA[2][0][1] );
	}

	BOOST_AUTO_TEST_CASE(multi_rotate_4d) {
		using std::get;  // workaround use of function template name with no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		multi::array<double, 4> original({14, 14, 7, 4});

		auto&& unrotd = original.unrotated();
		BOOST_TEST( get<0>(unrotd.sizes()) ==  4 );
		BOOST_TEST( get<1>(unrotd.sizes()) == 14 );
		BOOST_TEST( get<2>(unrotd.sizes()) == 14 );
		BOOST_TEST( get<3>(unrotd.sizes()) ==  7 );

		BOOST_TEST(( unrotd.sizes() == decltype(unrotd.sizes()){4, 14, 14, 7} ));
		BOOST_TEST( &original[0][1][2][3] == &unrotd[3][0][1][2] );

		auto&& unrotd2 = original.unrotated().unrotated();
		BOOST_TEST(( sizes(unrotd2) == decltype(sizes(unrotd2)){7, 4, 14, 14} ));
		BOOST_TEST( &original[0][1][2][3] == &unrotd2[2][3][0][1] );
	}

	BOOST_AUTO_TEST_CASE(multi_rotate_4d_op) {
		multi::array<double, 4> original({14, 14, 7, 4});

		auto&& unrotd = (original.unrotated());
		BOOST_TEST(( sizes(unrotd) == decltype(sizes(unrotd)){4, 14, 14, 7} ));
		BOOST_TEST( &original[0][1][2][3] == &unrotd[3][0][1][2] );

		auto&& unrotd2 = (original.unrotated().unrotated());
		BOOST_TEST(( sizes(unrotd2) == decltype(sizes(unrotd2)){7, 4, 14, 14} ));
		BOOST_TEST( &original[0][1][2][3] == &unrotd2[2][3][0][1] );
	}

	BOOST_AUTO_TEST_CASE(multi_rotate_part1) {
		// clang-format off
		std::array<std::array<int, 5>, 4> stdarr = {{
			{{ 0,  1,  2,  3,  4}},
			{{ 5,  6,  7,  8,  9}},
			{{10, 11, 12, 13, 14}},
			{{15, 16, 17, 18, 19}},
		}};
		// clang-format on

		std::array<std::array<int, 5>, 4> stdarr2 = {};

		multi::array_ref<int, 2> arr(&stdarr[0][0], {4, 5});    // NOLINT(readability-container-data-pointer) test access
		multi::array_ref<int, 2> arr2(&stdarr2[0][0], {4, 5});  // NOLINT(readability-container-data-pointer) test access

		arr2.rotated() = arr.rotated();

		BOOST_TEST( arr2[1][1] ==  6 );
		BOOST_TEST( arr2[2][1] == 11 );
		BOOST_TEST( arr2[1][2] ==  7 );

		BOOST_TEST( arr2.rotated()       == arr.rotated() );
		BOOST_TEST( arr2.rotated()[2][1] == 7             );
	}

	BOOST_AUTO_TEST_CASE(multi_rotate) {
		{
			multi::array<int, 2> arr = {
				{00, 01},
				{10, 11},
			};
			BOOST_TEST(       arr[1][0] == 10 );
			BOOST_TEST( (arr.rotated())[0][1] == 10 );
			BOOST_TEST( &     arr[1][0] == &(arr.rotated() )[0][1] );

			BOOST_TEST( arr.transposed()[0][1] == 10 );
			BOOST_TEST( arr.transposed()[0][1] == 10 );
			BOOST_TEST( (~arr)[0][1] == 10 );
			BOOST_TEST( &arr[1][0] == &arr.transposed()[0][1] );

			(arr.rotated())[0][1] = 100;
			BOOST_TEST( arr[1][0] == 100 );
		}
		{
			multi::array<double, 3> arr({11, 13, 17});
			BOOST_TEST( & arr[3][5][7] == &   arr   .transposed()[5][3][7] );
			BOOST_TEST( & arr[3][5][7] == &   arr   .transposed()[5][3][7] );
			BOOST_TEST( & arr[3][5][7] == & (~arr)            [5][3][7] );
			BOOST_TEST( & arr[3][5][7] == &   arr[3].transposed()[7][5] );
			BOOST_TEST( & arr[3][5][7] == & (~arr[3])            [7][5] );

			BOOST_TEST( & arr[3][5] == & (~arr)[5][3] );

			// BOOST_TEST( & ~~arr          == & arr      );
			// BOOST_TEST( &  (arr.rotated().rotated().rotated() )     == & arr       );
			// BOOST_TEST( &   arr()          == & (arr.rotated().rotated().rotated() ) );
			// BOOST_TEST( &  (arr.rotated() )     != & arr      );
			// BOOST_TEST( &  (arr.unrotated().rotated()) == & arr      );

			std::iota(arr.elements().begin(), arr.elements().end(), 0.1);
			BOOST_TEST( ~~arr == arr );
			BOOST_TEST( arr.unrotated().rotated() == arr );
		}
		{
			multi::array<int, 2> const arr = {
				{00, 01},
				{10, 11},
			};
			BOOST_TEST(   arr.rotated() [0][1] == 10 );
			BOOST_TEST( &(arr.rotated())[1][0] == &arr[0][1] );
			BOOST_TEST( &(~arr)[1][0] == &arr[0][1] );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_transposed) {
		multi::array<int, 2> const arr0 = {
			{ 9, 24, 30, 9},
			{ 4, 10, 12, 7},
			{14, 16, 36, 1},
		};
		multi::array<int, 2> const arr1 = arr0.transposed();
		multi::array<int, 2> const arr2 = ~arr0;
		BOOST_TEST( arr1 == arr2 );
	}

	BOOST_AUTO_TEST_CASE(miguel) {
		multi::array<double, 2> G2D({41, 35});
		auto const&             G3D = G2D.rotated().partitioned(7).sliced(0, 3).unrotated();

		BOOST_TEST( &G3D[0][0][0] == &G2D[0][0] );
	}
	{
		multi::array<int, 2> const arr = {
			{00, 01},
			{10, 11},
		};
		BOOST_TEST(   arr.rotated() [0][1] == 10 );
		BOOST_TEST( &(arr.rotated())[1][0] == &arr[0][1] );
		BOOST_TEST( &(~arr)[1][0] == &arr[0][1] );
	}

	BOOST_AUTO_TEST_CASE(multi_transposed) {
		multi::array<int, 2> const arr0 = {
			{ 9, 24, 30, 9},
			{ 4, 10, 12, 7},
			{14, 16, 36, 1},
		};
		multi::array<int, 2> const arr1{arr0.transposed()};
		multi::array<int, 2> const arr2{~arr0};
		BOOST_TEST( arr1 == arr2 );
	}

	BOOST_AUTO_TEST_CASE(miguel) {
		multi::array<double, 2> G2D({41, 35});
		auto const&             G3D = G2D.rotated().partitioned(7).sliced(0, 3).unrotated();

		BOOST_TEST( &G3D[0][0][0] == &G2D[0][0] );
	}

#if(__cplusplus >= 202002L)
	#if defined(__cpp_lib_ranges_repeat) && (__cpp_lib_ranges_repeat >= 202207L)
	#if !defined(__GNUC__) || (__GNUC__ < 14)
	BOOST_AUTO_TEST_CASE(matlab_meshgrid) {
		auto const x = multi::array{1, 2, 3};
		auto const y = multi::array{1, 2, 3, 4, 5};

		auto const& [X, Y] = meshgrid(x, y);

		auto const [X_copy, Y_copy] = meshgrid_copy(x, y);

		for(auto i : x.extension()) {
			for(auto j : y.extension()) {
				BOOST_TEST( X[i][j] == X_copy[i][j] );
				BOOST_TEST( Y[i][j] == Y_copy[i][j] );
			}
		}
	}
	#endif
	#endif
#endif
	return boost::report_errors();
}

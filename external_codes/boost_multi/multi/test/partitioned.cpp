// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, apply, subarray, operator==

#include <algorithm>  // for is_sorted
#include <array>      // for array
#include <cstddef>    // for ptrdiff_t
// #include <iostream>
#include <iterator>     // for size
#include <string>       // for operator""s, string, string_lite...
#include <tuple>        // for apply  // IWYU pragma: keep
#include <type_traits>  // for declval, decay_t, decay, decay<>...
#include <utility>      // for move

namespace multi = boost::multi;

template<class Ref> class propagate_const;

template<class T> class propagate_const<T&> {
	T& r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

 public:
	explicit propagate_const(T& other) : r_{other} {}
	propagate_const(propagate_const const&) = delete;
	propagate_const(propagate_const&&)      = delete;

	auto operator=(propagate_const const&) -> propagate_const&     = default;
	auto operator=(propagate_const&&) noexcept -> propagate_const& = default;

	auto operator=(T const& other) -> propagate_const& {
		r_ = other;
		return *this;
	}

	~propagate_const() noexcept = default;

	explicit operator T const&() const noexcept { return r_; }
	explicit operator T&() noexcept { return r_; }
};

template<class T> class propagate_const<T const&> {
	T const& r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

 public:
	explicit propagate_const(T const& other) : r_{other} {}
	auto     operator=(T const& other) -> propagate_const& = delete;
	explicit operator T const&() const noexcept { return r_; }
};

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(array_partitioned_1d) {
		multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50};

		auto&& A2_ref = A1.partitioned(2);

		static_assert(std::decay<decltype(A2_ref)>::type::rank{} == decltype(A1)::rank{} + 1);
		static_assert(std::decay_t<decltype(A2_ref)>::rank_v == decltype(A1)::rank_v + 1);

		BOOST_TEST( size(A2_ref   ) == 2 );
		BOOST_TEST( size(A2_ref[0]) == 3 );

		BOOST_TEST( &A2_ref[1][0] == &A1[3] );

		BOOST_TEST(( A2_ref == multi::array<double, 2>{ {0, 10, 20}, {30, 40, 50} } ));
	}

	BOOST_AUTO_TEST_CASE(array_partitioned_2d) {
		multi::array<int, 2> A2 = {
			{ 00,  10,  20,  30,  40,  50},
			{ 60,  70,  80,  90, 100, 110},

			{120, 130, 140, 150, 160, 170},
			{180, 190, 200, 210, 220, 230},
		};

		BOOST_TEST((
		A2.partitioned(2) == multi::array<int, 3>{
			{
				{ 00,  10,  20,  30,  40,  50},
				{ 60,  70,  80,  90, 100, 110},
			},
			{
				{120, 130, 140, 150, 160, 170},
				{180, 190, 200, 210, 220, 230},
			},
		}
	));

		auto&& A3_ref = A2.partitioned(2);

		static_assert(std::decay_t<decltype(A3_ref)>::rank{} == decltype(A2)::rank{} + 1);
		static_assert(std::decay_t<decltype(A3_ref)>::rank_v == decltype(A2)::rank_v + 1);

		BOOST_TEST( num_elements(A3_ref) == num_elements(A2) );
		BOOST_TEST( size(A3_ref) == 2 );
		BOOST_TEST( size(A3_ref[0]) == 2 );
		BOOST_TEST( size(A3_ref[0][0]) == 6 );
		BOOST_TEST( &A3_ref[1][1][0] == &A2[3][0] );

		A3_ref[0][0][0] = 99;
	}

	BOOST_AUTO_TEST_CASE(array_partitioned) {
		using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s
		using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		// NOLINTBEGIN(misc-include-cleaner) bug in clang-tidy 18
		multi::array<std::string, 2> A2 = {
			{"s0P0"s, "s1P0"s},
			{"s0P1"s, "s1P1"s},
			{"s0P2"s, "s1P2"s},
			{"s0P3"s, "s1P3"s},
			{"s0P4"s, "s1P4"s},
			{"s0P5"s, "s1P5"s},
		};
		// NOLINTEND(misc-include-cleaner)

		BOOST_TEST(  size(A2) == 6 );

		BOOST_TEST( get<0>(A2.sizes()) == 6 );
		BOOST_TEST( get<1>(A2.sizes()) == 2 );

		BOOST_TEST(( A2.sizes() == decltype(A2.sizes()){6, 2} ));

		BOOST_TEST( get<0>(A2.sizes()) == 6 );
		BOOST_TEST( get<1>(A2.sizes()) == 2 );

		BOOST_TEST( size(A2.partitioned(3)) == 3 );

		static_assert(decltype(A2.partitioned(3))::rank{} == 3);
		static_assert(decltype(A2.partitioned(3))::rank::value == 3);
		static_assert(decltype(A2.partitioned(3))::rank_v == 3);

		BOOST_TEST(( sizes(A2.partitioned(3)) == decltype(sizes(A2.partitioned(3))){3, 2, 2} ));

		BOOST_TEST( get<0>(sizes(A2.partitioned(3))) == 3 );
		BOOST_TEST( get<1>(sizes(A2.partitioned(3))) == 2 );
		BOOST_TEST( get<2>(sizes(A2.partitioned(3))) == 2 );

		BOOST_TEST( size(A2.partitioned(1)) == 1 );

		static_assert(decltype(A2.partitioned(1))::rank{} == 3);
		static_assert(decltype(A2.partitioned(1))::rank::value == 3);
		static_assert(decltype(A2.partitioned(1))::rank_v == 3);

		BOOST_TEST( &A2.partitioned(1).rotated()[3][1][0] == &A2[3][1] );
	}

	BOOST_AUTO_TEST_CASE(array_encoded_subarray) {
		using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		// arr[walker][encoded_property]  // 7 walkers
		multi::array<int, 2> arr = {
			{990, 990, 1000, 001,  10,  11,  20,  21, 990},
			{990, 990,  100, 101, 110, 111, 120, 121, 990},
			{990, 990,  200, 201, 210, 211, 220, 221, 990},
			{990, 990,  300, 301, 310, 311, 320, 321, 990},
			{990, 990,  400, 401, 410, 411, 420, 421, 990},
			{990, 990,  500, 501, 510, 511, 520, 521, 990},
			{990, 990,  600, 601, 610, 611, 620, 621, 990},
		};

		// multi::iextension const encoded_3x2_range = {2, 8};

		auto&& arrRPU = arr.rotated().sliced(2, 8).partitioned(3).unrotated();

		static_assert(decltype(+arrRPU)::rank::value == 3);
		static_assert(decltype(+arrRPU)::rank{} == 3);
		static_assert(decltype(+arrRPU)::rank_v == 3);

		BOOST_TEST( get<0>(arrRPU.sizes()) == 7 );
		BOOST_TEST( get<1>(arrRPU.sizes()) == 3 );
		BOOST_TEST( get<2>(arrRPU.sizes()) == 2 );

		BOOST_TEST(( arrRPU.sizes() == decltype(arrRPU.sizes()){7, 3, 2} ));
		BOOST_TEST(( sizes(arrRPU) == decltype(sizes(arrRPU)){7, 3, 2} ));
		BOOST_TEST( arrRPU[4].num_elements() == 3*2L );

		BOOST_TEST( &arrRPU[4][1][0] == &arr[4][4] );
		BOOST_TEST( arrRPU[4][1][0] == 410 );

		BOOST_TEST((
		arrRPU[4] == multi::array<double, 2>{
			{400, 401},
			{410, 411},
			{420, 421},
		}
	));

		arrRPU[4][1][0] = 11110;
		BOOST_TEST( arr[4][4] == 11110 );

		class walker_ref {
			using raw_source_reference = decltype(std::declval<multi::array<int, 2>&>()[0]);
			using internal_array_type  = decltype(std::declval<raw_source_reference>().sliced(2, 8).partitioned(3));

		 public:                                 // NOLINT(whitespace/indent) bug in cpplint
			propagate_const<int&> prop1;         // NOLINT(misc-non-private-member-variables-in-classes)
			propagate_const<int&> prop2;         // NOLINT(misc-non-private-member-variables-in-classes)
			internal_array_type   slater_array;  // NOLINT(misc-non-private-member-variables-in-classes)
			propagate_const<int&> prop3;         // NOLINT(misc-non-private-member-variables-in-classes)

			explicit walker_ref(raw_source_reference&& row)
			: prop1{row[0]}, prop2{row[1]}, slater_array{row.sliced(2, 8).partitioned(3)}, prop3{std::move(row)[8]} {}
		};

		auto&& wr = walker_ref(arr[5]);

		wr.prop1 = 88;

		BOOST_TEST( wr.slater_array[2][1] == 521 );

		// what( wr , wr.slater_array, wr.slater_array[2][1] );
		wr.slater_array[2][1] = 99990;
	}

	BOOST_AUTO_TEST_CASE(array_partitioned_add_to_last) {
		multi::array<double, 3> arr = {
			{
             {0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
             {6.0, 7.0, 8.0, 9.0, 10.0, 11.0},
             {12.0, 13.0, 14.0, 15.0, 16.0, 17.0},
             {18.0, 19.0, 20.0, 21.0, 22.0, 23.0},
			 },
			{
             {0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
             {6.0, 7.0, 8.0, 9.0, 10.0, 11.0},
             {12.0, 13.0, 14.0, 15.0, 16.0, 17.0},
             {18.0, 19.0, 20.0, 21.0, 22.0, 23.0},
			 },
		};

		auto strides = std::apply([](auto... strds) { return std::array<std::ptrdiff_t, sizeof...(strds)>{{strds...}}; }, arr.layout().strides());
		// auto strides = std::apply([](auto... strds) { return std::array<std::ptrdiff_t, sizeof...(strds)>{{strds...}}; }, arr.strides());

		BOOST_TEST( std::is_sorted(strides.rbegin(), strides.rend()) && arr.num_elements() == arr.nelems() );  // contiguous c-ordering

#ifndef _MSC_VER  // problem with MSVC 14.3 c++17
		auto&& A4 = arr.reinterpret_array_cast<double>(1);

		BOOST_TEST(( arr.extensions() == decltype(arr.extensions()){2, 4, 6} ));
		BOOST_TEST(( A4.extensions() == decltype(A4.extensions()){2, 4, 6, 1} ));

		//  BOOST_TEST( A4.is_flattable() );
		//  BOOST_TEST( A4.flatted().is_flattable() );

		BOOST_TEST( &A4[1][2][3][0] == &arr[1][2][3] );
#endif
	}

	BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_1D) {
		multi::array<double, 1> arr = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
		BOOST_TEST( size(arr.partitioned(3)) == 3 );
		BOOST_TEST(( arr.partitioned(3)[1] == multi::array<double, 1>{4.0, 5.0, 6.0, 7.0} ));
		BOOST_TEST( &arr.partitioned(3)[1][2] == &arr[6] );

		BOOST_TEST( size(arr.chunked(3)) == 4 );
		BOOST_TEST(( arr.chunked(3)[1] == multi::array<double, 1>({3.0, 4.0, 5.0}) ));
		BOOST_TEST( &arr.chunked(3)[1][2] == &arr[5] );
	}

	BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_2D) {
		multi::array<double, 2> arr({100, 53});
		BOOST_TEST( arr.partitioned(20).size() == 20 );
		BOOST_TEST( &arr.partitioned(20)[1][2] == &arr[7] );

		BOOST_TEST( size(arr.chunked(5)) == 20 );
		BOOST_TEST( &arr.chunked(5)[1][2] == &arr[7] );
	}

	BOOST_AUTO_TEST_CASE(chunked_subarrays) {
		using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		multi::array<int, 2> const arr = {
			{ 0,  1, /**/ 2,   3, /**/  4,  5},
			{ 6,  7, /**/ 8,   9, /**/ 10, 11},
			/*********************************/
			{12, 13, /**/ 14, 15, /**/ 16, 17},
			{18, 19, /**/ 20, 21, /**/ 22, 23},
			/*********************************/
			{24, 25, /**/ 26, 27, /**/ 28, 29},
			{30, 31, /**/ 32, 33, /**/ 34, 35},
			/*********************************/
			{36, 37, /**/ 38, 39, /**/ 40, 41},
			{42, 43, /**/ 44, 45, /**/ 46, 47}
		};

		BOOST_TEST( arr.dimensionality == 2 );
		BOOST_TEST( arr.size() == 8 );
		BOOST_TEST( get<1>(arr.sizes()) == 6 );

		BOOST_TEST( arr.chunked(2).dimensionality == 3 );
		BOOST_TEST( arr.chunked(2).size() == 4 );
		BOOST_TEST( get<1>(arr.chunked(2).sizes()) == 2 );
		BOOST_TEST( get<2>(arr.chunked(2).sizes()) == 6 );

		auto const&& block_arr = arr.chunked(2).rotated().rotated().chunked(2).transposed().rotated().transposed();
		BOOST_TEST( block_arr.dimensionality == 4 );

		BOOST_TEST( block_arr.size() == 4 );
		BOOST_TEST( get<1>(block_arr.sizes()) == 3 );
		BOOST_TEST( get<2>(block_arr.sizes()) == 2 );
		BOOST_TEST( get<3>(block_arr.sizes()) == 2 );

		BOOST_TEST((
			block_arr[2][1]
			==
				multi::array<int, 2>{
					{26, 27},
					{32, 33}
				}
		));
	}

	BOOST_AUTO_TEST_CASE(partitined_subarrays) {
		using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		multi::array<int, 2> const arr = {
			{ 0,  1, /**/ 2,   3, /**/  4,  5},
			{ 6,  7, /**/ 8,   9, /**/ 10, 11},
			/******************************************/
			{12, 13, /**/ 14, 15, /**/ 16, 17},
			{18, 19, /**/ 20, 21, /**/ 22, 23},
			/******************************************/
			{24, 25, /**/ 26, 27, /**/ 28, 29},
			{30, 31, /**/ 32, 33, /**/ 34, 35},
			/******************************************/
			{36, 37, /**/ 38, 39, /**/ 40, 41},
			{42, 43, /**/ 44, 45, /**/ 46, 47}
		};

		BOOST_TEST( arr.dimensionality == 2 );
		BOOST_TEST( arr.partitioned(4).dimensionality == 3 );

		BOOST_TEST( arr.partitioned(4).size() == 4 );
		BOOST_TEST( get<1>(arr.partitioned(4).sizes()) == 2 );
		BOOST_TEST( get<2>(arr.partitioned(4).sizes()) == 6 );

		auto const&& block_arr = arr.partitioned(4).rotated().rotated().partitioned(3).transposed().rotated().transposed();
		BOOST_TEST( block_arr.dimensionality == 4 );

		BOOST_TEST( block_arr.size() == 4 );
		BOOST_TEST( get<1>(block_arr.sizes()) == 3 );
		BOOST_TEST( get<2>(block_arr.sizes()) == 2 );
		BOOST_TEST( get<3>(block_arr.sizes()) == 2 );

		BOOST_TEST((
			block_arr[2][1]
			==
			multi::array<int, 2>{
				{26, 27},
				{32, 33}
			}
		));
	}

	BOOST_AUTO_TEST_CASE(tiled_1D) {
		multi::array<int, 1> const arr = {0, 1, 2, /**/ 3, 4, 5, /**/ 6, 7};

		BOOST_TEST( arr.size() == 8 );

		BOOST_TEST(( arr.tiled(3).quotient[0].size() == 3 ));
		BOOST_TEST(( arr.tiled(3).quotient[1].size() == 3 ));

		BOOST_TEST(( arr.tiled(3).remainder.size() == 2 ));

		auto [tiles, border] = arr.tiled(3);
		BOOST_TEST( tiles.size() == 2 );
		BOOST_TEST(( tiles[0].size() == 3 ));
		BOOST_TEST(( tiles[1].size() == 3 ));

		BOOST_TEST(( tiles[0] == multi::array<int, 1>{0, 1, 2} ));
		BOOST_TEST(( tiles[1] == multi::array<int, 1>{3, 4, 5} ));

		BOOST_TEST( border.size() == 2 );
		BOOST_TEST(( border == multi::array<int, 1>{6, 7} ));
	}

	return boost::report_errors();
}

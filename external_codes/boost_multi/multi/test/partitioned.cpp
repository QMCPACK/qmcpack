// Copyright 2018-2026 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, apply, subarray, operator==

#include <boost/core/lightweight_test.hpp>

#include <algorithm>    // for is_sorted
#include <array>        // for array
#include <cstddef>      // for ptrdiff_t
#include <string>       // for operator""s, string, string_lite...
#include <tuple>        // for apply  // IWYU pragma: keep
#include <type_traits>  // for declval, decay_t, decay, decay<>...
#include <utility>      // for move

#ifdef _MSC_VER
#pragma warning(disable : 4625)  // copy constructor was implicitly defined as deleted
#pragma warning(disable : 4626)  // assignment operator was implicitly defined as deleted
#pragma warning(disable : 5026)  // move constructor was implicitly defined as deleted
#endif

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

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// // BOOST_AUTO_TEST_CASE(halved_1d)
	// {
	// 	multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50};

	// 	auto&& A2_ref = A1.halved();

	// 	// static_assert(std::decay_t<decltype(A2_ref)>::rank{} == decltype(A1)::rank{} + 1);
	// 	static_assert(std::decay_t<decltype(A2_ref)>::rank_v == decltype(A1)::rank_v + 1);

	// 	BOOST_TEST( A2_ref.size() == 2 );
	// 	BOOST_TEST( A2_ref[0].size() == 3 );

	// 	BOOST_TEST( &A2_ref[1][0] == &A1[3] );

	// 	BOOST_TEST(( A2_ref == multi::array<double, 2>{ {0, 10, 20}, {30, 40, 50} } ));
	// }

	// // BOOST_AUTO_TEST_CASE(halved_2d)
	// {
	// 	multi::array<int, 2> const A2 = {
	// 		{ 00,  10,  20,  30,  40,  50},
	// 		{ 60,  70,  80,  90, 100, 110},

	// 		{120, 130, 140, 150, 160, 170},
	// 		{180, 190, 200, 210, 220, 230},
	// 	};

	// 	BOOST_TEST((
	// 		A2.halved() == multi::array<int, 3>{
	// 			{
	// 				{ 00,  10,  20,  30,  40,  50},
	// 				{ 60,  70,  80,  90, 100, 110},
	// 			},
	// 			{
	// 				{120, 130, 140, 150, 160, 170},
	// 				{180, 190, 200, 210, 220, 230},
	// 			},
	// 		}
	// 	));
	// }

	// // BOOST_AUTO_TEST_CASE(halved_ref_2d)
	// {
	// 	std::vector<int> buff({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}, {});
	// 	auto const&      arr = multi::array_ref<int, 2>({4, 4}, buff.data());

	// 	BOOST_TEST(( arr == multi::array<int, 2>{
	// 		{ 1,  2,  3,  4},
	// 		{ 5,  6,  7,  8},
	// 		{ 9, 10, 11, 12},
	// 		{13, 14, 15, 16}
	// 	}));

	// 	BOOST_TEST((  arr.halved().extents() == multi::extents_t{2, 2, 4} ));
	// 	BOOST_TEST(( arr.halved() == multi::array<int, 3>{
	// 		{
	// 			{ 1,  2,  3,  4},
	// 			{ 5,  6,  7,  8}
	// 		},
	// 		{
	// 			{ 9, 10, 11, 12},
	// 			{13, 14, 15, 16}
	// 		}
	// 	}));

	// 	BOOST_TEST(( arr.halved().flatted() == multi::array<int, 2>{
	// 		{ 1,  2,  3,  4},
	// 		{ 5,  6,  7,  8},
	// 		{ 9, 10, 11, 12},
	// 		{13, 14, 15, 16}
	// 	}));

	// 	BOOST_TEST(( arr.halved().transposed() == multi::array<int, 3>{
	// 		{
	// 			{ 1,  2,  3,  4},
	// 			{ 9, 10, 11, 12},
	// 		},
	// 		{
	// 			{ 5,  6,  7,  8},
	// 			{13, 14, 15, 16}
	// 		}
	// 	}));

	// 	BOOST_TEST(( arr.halved().rotated() == multi::array<int, 3>{
	// 		{
	// 			{1, 9},
	// 			{2, 10},
	// 			{3, 11},
	// 			{4, 12}
	// 		},
	// 		{
	// 			{ 5, 13},
	// 			{ 6, 14},
	// 			{ 7, 15},
	// 			{ 8, 16}
	// 		}
	// 	}));

	// 	BOOST_TEST(( arr.halved().rotated().rotated() == multi::array<int, 3>{
	// 		{
	// 			{1, 5},
	// 			{9, 13}
	// 		},
	// 		{
	// 			{2, 6},
	// 			{10, 14}
	// 		},
	// 		{
	// 			{3, 7},
	// 			{11, 15}
	// 		},
	// 		{
	// 			{4, 8},
	// 			{12, 16}
	// 		}
	// 	}));

	// 	BOOST_TEST(( arr.halved().rotated().rotated().halved() == multi::array<int, 4>{
	// 		{
	// 			{
	// 				{ 1,  5},
	// 				{ 9, 13}
	// 			},
	// 			{
	// 				{ 2,  6},
	// 				{10, 14}
	// 			}
	// 		},
	// 		{
	// 			{
	// 				{ 3,  7},
	// 				{11, 15}
	// 			},
	// 			{
	// 				{ 4,  8},
	// 				{12, 16}
	// 			}
	// 		}
	// 	}));

	// 	BOOST_TEST(( arr.halved().rotated().rotated().halved().unrotated() == multi::array<int, 4>{
	// 		{
	// 			{
	// 				{ 1,  9},
	// 				{ 2,  10}
	// 			},
	// 			{
	// 				{ 3, 11},
	// 				{ 4, 12}
	// 			}
	// 		},
	// 		{
	// 			{
	// 				{ 5, 13},
	// 				{ 6, 14}
	// 			},
	// 			{
	// 				{ 7, 15},
	// 				{ 8, 16}
	// 			}
	// 		}
	// 	}));

	// 	BOOST_TEST( arr.halved().rotated().rotated().halved().unrotated().unrotated()[0][0][0][0] == 1 );
	// 	BOOST_TEST( arr.halved().rotated().rotated().halved().unrotated().unrotated()[0][0][0][1] == 2 );
	// 	BOOST_TEST( arr.halved().rotated().rotated().halved().unrotated().unrotated()[0][0][1][0] == 3 );
	// 	BOOST_TEST( arr.halved().rotated().rotated().halved().unrotated().unrotated()[0][0][1][1] == 4 );
	// 	BOOST_TEST( arr.halved().rotated().rotated().halved().unrotated().unrotated()[0][1][0][0] == 5 );
	// 	BOOST_TEST( arr.halved().rotated().rotated().halved().unrotated().unrotated()[0][1][0][1] == 6 );

	// 	BOOST_TEST(( arr.halved().rotated().rotated().halved().unrotated().unrotated() == multi::array<int, 4>{
	// 		{
	// 			{
	// 				{ 1,  2},
	// 				{ 3,  4}
	// 			},
	// 			{
	// 				{ 5, 6},
	// 				{ 7, 8}
	// 			}
	// 		},
	// 		{
	// 			{
	// 				{ 9, 10},
	// 				{11, 12}
	// 			},
	// 			{
	// 				{13, 14},
	// 				{15, 16}
	// 			}
	// 		}
	// 	}));

	// 	BOOST_TEST(( arr.halved().rotated().rotated().halved().unrotated().unrotated().unrotated() == multi::array<int, 4>{
	// 		{
	// 			{
	// 				{ 1,  3},
	// 				{ 5,  7}
	// 			},
	// 			{
	// 				{ 9, 11},
	// 				{13, 15}
	// 			}
	// 		},
	// 		{
	// 			{
	// 				{ 2, 4},
	// 				{ 6, 8}
	// 			},
	// 			{
	// 				{10, 12},
	// 				{14, 16}
	// 			}
	// 		}
	// 	}));

	// 	BOOST_TEST(( arr.halved().rotated().rotated().halved().unrotated().unrotated().unrotated().unrotated() == multi::array<int, 4>{
	// 		{
	// 			{
	// 				{ 1,  5},
	// 				{ 9, 13}
	// 			},
	// 			{
	// 				{ 2, 6},
	// 				{10, 14}
	// 			}
	// 		},
	// 		{
	// 			{
	// 				{ 3, 7},
	// 				{11, 15}
	// 			},
	// 			{
	// 				{ 4,  8},
	// 				{12, 16}
	// 			}
	// 		}
	// 	}));
	// }

	// BOOST_AUTO_TEST_CASE(array_partitioned_1d)
	{
		multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50};

		auto&& A2_ref = A1.partitioned(2);

		auto [is, js] = A2_ref.extents();
		BOOST_TEST( is.size() == 2 );
		BOOST_TEST( js.size() == 3 );

		// static_assert(std::decay<decltype(A2_ref)>::type::rank{} == decltype(A1)::rank{} + 1);
		static_assert(std::decay_t<decltype(A2_ref)>::rank_v == decltype(A1)::rank_v + 1);

		BOOST_TEST( A2_ref.size() == 2 );
		BOOST_TEST( A2_ref[0].size() == 3 );

		BOOST_TEST( &A2_ref[1][0] == &A1[3] );

		BOOST_TEST(( A2_ref == multi::array<double, 2>{ {0, 10, 20}, {30, 40, 50} } ));
	}

	// BOOST_AUTO_TEST_CASE(array_partitioned_2d)
	{
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

		static_assert(std::decay_t<decltype(A3_ref)>::dimensionality == decltype(A2)::dimensionality + 1);

		BOOST_TEST( A3_ref.num_elements() == A2.num_elements() );
		BOOST_TEST( A3_ref.size() == 2 );
		BOOST_TEST( A3_ref[0].size() == 2 );
		BOOST_TEST( A3_ref[0][0].size() == 6 );
		BOOST_TEST( &A3_ref[1][1][0] == &A2[3][0] );

		A3_ref[0][0][0] = 99;
		BOOST_TEST( A3_ref[0][0][0] == 99 );  // cppcheck-suppress knownConditionTrueFalse ;
	}

	// BOOST_AUTO_TEST_CASE(array_partitioned)
	{
		using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s
		using std::get;                        // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

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

		BOOST_TEST( A2.partitioned(3).size() == 3 );

		static_assert(decltype(A2.partitioned(3))::dimensionality == 3);
		// static_assert(decltype(A2.partitioned(3))::rank{} == 3);
		// static_assert(decltype(A2.partitioned(3))::rank::value == 3);
		// static_assert(decltype(A2.partitioned(3))::rank_v == 3);

		BOOST_TEST(( sizes(A2.partitioned(3)) == decltype(sizes(A2.partitioned(3))){3, 2, 2} ));

		BOOST_TEST( get<0>(sizes(A2.partitioned(3))) == 3 );
		BOOST_TEST( get<1>(sizes(A2.partitioned(3))) == 2 );
		BOOST_TEST( get<2>(sizes(A2.partitioned(3))) == 2 );

		BOOST_TEST( A2.partitioned(1).size() == 1 );

		static_assert(decltype(A2.partitioned(1))::dimensionality == 3);

		BOOST_TEST( &A2.partitioned(1).rotated()[3][1][0] == &A2[3][1] );
	}

	// BOOST_AUTO_TEST_CASE(array_encoded_subarray)
	{
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

		static_assert(decltype(+arrRPU)::dimensionality == 3);

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

		wr.slater_array[2][1] = 99990;
		BOOST_TEST( wr.slater_array[2][1] == 99990 );
	}

	// BOOST_AUTO_TEST_CASE(array_partitioned_add_to_last)
	{
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

		// NOLINTNEXTLINE(modernize-use-ranges) for C++20
		BOOST_TEST( std::is_sorted(strides.rbegin(), strides.rend()) && arr.num_elements() == arr.nelems() );  // contiguous c-ordering

		// #ifndef _MSC_VER  // problem with MSVC 14.3 c++17
		auto&& A4 = arr.reinterpret_array_cast<double>(1);

		BOOST_TEST(( arr.extents() == decltype(arr.extents()){2, 4, 6} ));
		BOOST_TEST(( A4.extents() == decltype(A4.extents()){2, 4, 6, 1} ));

		//  BOOST_TEST( A4.is_flattable() );
		//  BOOST_TEST( A4.flatted().is_flattable() );

		BOOST_TEST( &A4[1][2][3][0] == &arr[1][2][3] );
		// #endif
	}

	// BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_1D)
	{
		multi::array<double, 1> arr = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
		BOOST_TEST( arr.partitioned(3).size() == 3 );
		BOOST_TEST(( arr.partitioned(3)[1] == multi::array<double, 1>{4.0, 5.0, 6.0, 7.0} ));
		BOOST_TEST( &arr.partitioned(3)[1][2] == &arr[6] );

		BOOST_TEST(  arr.chunked(3).size() == 4 );
		BOOST_TEST(( arr.chunked(3)[1] == multi::array<double, 1>({3.0, 4.0, 5.0}) ));
		BOOST_TEST( &arr.chunked(3)[1][2] == &arr[5] );
	}

	// BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_2D)
	{
		multi::array<double, 2> arr({100, 53});
		BOOST_TEST( arr.partitioned(20).size() == 20 );
		BOOST_TEST( &arr.partitioned(20)[1][2] == &arr[7] );

		BOOST_TEST( arr.chunked(5).size() == 20 );
		BOOST_TEST( &arr.chunked(5)[1][2] == &arr[7] );
	}

	// BOOST_AUTO_TEST_CASE(chunked_subarrays)
	{
		using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		multi::array<int, 2> const arr = {
			{ 0,  1,  /**/ 2,  3,  /**/ 4,  5},
			{ 6,  7,  /**/ 8,  9, /**/ 10, 11},
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

	// BOOST_AUTO_TEST_CASE(partitined_subarrays)
	{
		multi::array<int, 2> const arr = {
			{ 0,  1,  /**/ 2,  3,  /**/ 4,  5},
			{ 6,  7,  /**/ 8,  9, /**/ 10, 11},
			/******************************************/
			{12, 13, /**/ 14, 15, /**/ 16, 17},
			{18, 19, /**/ 20, 21, /**/ 22, 23},
			/******************************************/
			{24, 25, /**/ 26, 27, /**/ 28, 29},
			{30, 31, /**/ 32, 33, /**/ 34, 35},
			/******************************************/
			{36, 37, /**/ 38, 39, /**/ 40, 41},
			{42, 43, /**/ 44, 45, /**/ 46, 47},
		};

		BOOST_TEST( arr.dimensionality == 2 );
		BOOST_TEST( arr.partitioned(4).dimensionality == 3 );

		BOOST_TEST( arr.partitioned(4).size() == 4 );

		using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

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

	// tiled 1D with 3 elements per tile
	{
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

	// tiled 1D with 5 elements per tile
	{
		multi::array<int, 1> const arr = {0, 1, 2, 3, 4, /**/ 5, 6, 7, 8, 9, /**/ 10, 11, 12, 13};

		BOOST_TEST( arr.size() == 14 );

		BOOST_TEST(( arr.tiled(5).quotient[0].size() == 5 ));
		BOOST_TEST(( arr.tiled(5).quotient[1].size() == 5 ));

		BOOST_TEST(( arr.tiled(5).remainder.size() == 4 ));

		auto [tiles, border] = arr.tiled(5);
		BOOST_TEST( tiles.size() == 2 );
		BOOST_TEST(( tiles[0].size() == 5 ));
		BOOST_TEST(( tiles[1].size() == 5 ));

		BOOST_TEST(( tiles[0] == multi::array<int, 1>{0, 1, 2, 3, 4} ));
		BOOST_TEST(( tiles[1] == multi::array<int, 1>{5, 6, 7, 8, 9} ));

		BOOST_TEST( border.size() == 4 );
		BOOST_TEST(( border == multi::array<int, 1>{10, 11, 12, 13} ));
	}

	// put subarray in std::array
	{
		multi::array<int, 1> arr = {1, 2, 3, /**/ 4, 5};
		BOOST_TEST( arr.size() == 5 );

		// std::array<multi::array<int, 1>::subarray, 2> arr2 = {{ arr({0, 3}), arr({3, 5}) }};
		// auto arr2 = std::array<multi::array<int, 1>::subarray, 2>{{ arr({0, 3}), arr({3, 5}) }};
		// auto&& arr2 = std::array<multi::array<int, 1>::subarray, 2>{{ arr({0, 3}), arr({3, 5}) }};

		// CTAD fails:
		// auto&& arr2 = std::array{{ arr({0, 3}), arr({3, 5}) }};

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
#endif
		auto&& arr2 = std::array{arr({0, 3}), arr({3, 5})};
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

		// arr2 is not copyable, good
		// auto arr3 = arr2;
		[[maybe_unused]] auto&& arr3 = arr2;

		BOOST_TEST( arr2[0].size() == 3 );
		BOOST_TEST( arr2[1].size() == 2 );

		BOOST_TEST( &arr2[0][1] == &arr[1] );
		BOOST_TEST( &arr2[1][1] == &arr[4] );

		arr2[0][1] = 99;
		BOOST_TEST( arr[1] == 99 );
	}

	// put subarrays in tuple
	{
		multi::array<int, 1> arr = {1, 2, 3, /**/ 4, 5};
		BOOST_TEST( arr.size() == 5 );

		// asan complains in these cases:
		// std::tuple<multi::array<int, 1>::subarray&&, multi::array<int, 1>::subarray&&> tup(arr({0, 3}), arr({3, 5}));
		// auto tup = std::tuple<multi::array<int, 1>::subarray&&, multi::array<int, 1>::subarray&&>(arr({0, 3}), arr({3, 5}));

		// invalid syntax:
		// std::tuple<multi::array<int, 1>::subarray&&, multi::array<int, 1>::subarray&&>&& tup(arr({0, 3}), arr({3, 5}));

		// asan complains:
		// auto&& tup = std::tuple<multi::array<int, 1>::subarray&&, multi::array<int, 1>::subarray&&>(arr({0, 3}), arr({3, 5}));

		// asan complains:
		// auto&& tup = std::forward_as_tuple(arr({0, 3}), arr({3, 5}));
		// std::tuple<multi::array<int, 1>::subarray, multi::array<int, 1>::subarray> tup(arr({0, 3}), arr({3, 5}));
		// auto tup = std::tuple<multi::array<int, 1>::subarray, multi::array<int, 1>::subarray>(arr({0, 3}), arr({3, 5}));
		// auto tup = std::make_tuple(arr({0, 3}), arr({3, 5}));
		// auto tup = std::tuple(arr({0, 3}), arr({3, 5}));
		// auto&& tup = std::tuple(arr({0, 3}), arr({3, 5}));

		auto&& tup = std::tuple{arr({0, 3}), arr({3, 5})};

		// tup is not copyable, good
		// auto tup2 = tup;
		[[maybe_unused]] auto&& tup2 = tup;

		using std::get;  // for C++17
		BOOST_TEST( get<0>(tup).size() == 3 );
		BOOST_TEST( get<1>(tup).size() == 2 );

		BOOST_TEST( &get<0>(tup)[1] == &arr[1] );
		BOOST_TEST( &get<1>(tup)[1] == &arr[4] );

		get<0>(tup)[1] = 99;
		BOOST_TEST( arr[1] == 99 );
	}

	// put 2D subarrays in 2D tuple
	{
		multi::array<int, 2> arr = {
			{1, 2, 3, /**/ 4,  5},
			{6, 7, 8, /**/ 9, 10},
			/*******************/
			{1, 2, 3, /**/ 4,  5},
			{6, 7, 8, /**/ 9, 10},
		};

		BOOST_TEST( arr.elements().size() == 20 );

		auto&& tup = std::tuple{
			arr({0, 2}, {0, 3}),
			arr({0, 2}, {3, 5}),
			arr({2, 4}, {0, 3}),
			arr({2, 4}, {3, 5}),
		};

		using std::get;  // for C++17
		BOOST_TEST( &get<0>(tup)[1][1] == &arr[1][1] );
		BOOST_TEST( &get<2>(tup)[0][2] == &arr[2][2] );

		auto&& tup2d = std::tuple{
			std::tuple{arr({0, 2}, {0, 3}), arr({0, 2}, {3, 5})},
			std::tuple{arr({2, 4}, {0, 3}), arr({2, 4}, {3, 5})},
		};

		BOOST_TEST( &get<0>(get<0>(tup2d))[1][1] == &arr[1][1] );
		BOOST_TEST( &get<0>(get<1>(tup2d))[0][2] == &arr[2][2] );
	}

	// put 2D subarrays in 2D tuple
	{
		multi::array<int, 2> arr = {
			{1, 2, 3, /**/ 4,  5},
			{6, 7, 8, /**/ 9, 10},
			/*******************/
			{1, 2, 3, /**/ 4,  5},
			{6, 7, 8, /**/ 9, 10},
		};

		BOOST_TEST( arr.elements().size() == 20 );

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
#endif
		auto&& tup = std::array{
			arr({0, 2}, {0, 3}),
			arr({0, 2}, {3, 5}),
			arr({2, 4}, {0, 3}),
			arr({2, 4}, {3, 5}),
		};

		using std::get;  // for C++17
		BOOST_TEST( &get<0>(tup)[1][1] == &arr[1][1] );
		BOOST_TEST( &get<2>(tup)[0][2] == &arr[2][2] );

#if defined(__clang__) || (!defined(__GNUC__) || (__GNUC__ > 9))  // gcc 9 gets confused with CTAD
		auto&& tup2d = std::array{
			std::array{arr({0, 2}, {0, 3}), arr({0, 2}, {3, 5})},
			std::array{arr({2, 4}, {0, 3}), arr({2, 4}, {3, 5})},
		};
		BOOST_TEST( &tup2d[0][0] [1][1] == &arr[1][1] );
		BOOST_TEST( &tup2d[1][0] [0][2] == &arr[2][2] );
#endif
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
	}
	// {
	// 	multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50};

	// 	auto&& [left, right] = A1.splitted();

	// 	BOOST_TEST( left.size() == 3 );
	// 	BOOST_TEST( right.size() == 3 );

	// 	BOOST_TEST( A1[1] == 10 );
	// }
	// {
	// 	multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50};

	// 	auto&& [left, right] = std::move(A1).splitted();

	// 	BOOST_TEST( left.size() == 3 );
	// 	BOOST_TEST( right.size() == 3 );

	// 	// BOOST_TEST( A1[1] == 10 ); use after move detected by clang-tidy bugprone-use-after-move
	// }
	// {
	// 	multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50};

	// 	auto&& [left, right] = A1.split();

	// 	BOOST_TEST( left.size() == 3 );
	// 	BOOST_TEST( right.size() == 3 );

	// 	// BOOST_TEST( A1[1] == 10 ); use after move detected by clang-tidy bugprone-use-after-move
	// }
	// {
	// 	multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50, 60};

	// 	auto&& [left, right] = A1.splitted();

	// 	BOOST_TEST( left.size() == 3 );
	// 	BOOST_TEST( right.size() == 4 );
	// }
	// {
	// 	multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50, 60, 70};

	// 	auto const& A2 = A1.strided(2);
	// 	BOOST_TEST( A2.size() == 4 );

	// 	auto&& [left, right] = A2.splitted();

	// 	BOOST_TEST( left.size() == 2 );
	// 	BOOST_TEST( right.size() == 2 );
	// }
	// {
	// 	multi::array<int, 1> A1 = {0, 10, 20, 30, 40, 50, 60, 70};

	// 	auto const& A2 = A1.strided(2);
	// 	BOOST_TEST( A2.size() == 4 );

	// 	auto&& [left, right] = A2.splitted();

	// 	BOOST_TEST( left.size() == 2 );
	// 	BOOST_TEST( right.size() == 2 );
	// }

	// -Wconsumed probe: split() is annotated [[clang::set_typestate(consumed)]],
	// so using `A` after the call must trigger a clang -Wconsumed warning.
	// Behavior gate:
	//   - default                  : noisy warning (informational)
	//   - -DMULTI_EXPECT_CONSUMED_WARNING with -Werror=consumed : turns into a
	//     hard error if the warning ever stops firing (regression detector)
	// {
	// 	multi::array<int, 1> arr = {0, 10, 20, 30};

	// 	auto&& [left, right] = arr.split();  // A is now in "consumed" typestate
	// 	BOOST_TEST( left.size() == 2 );
	// 	BOOST_TEST( right.size() == 2 );

	// 	// any use of A here should fire -Wconsumed; we exercise a few:
	// 	BOOST_TEST( arr.size() == 4 );  // expected: warning: invalid invocation of method 'size' on object 'A' while it is in the 'consumed' state
	// 	BOOST_TEST( arr[0] == 0 );      // expected: warning: invalid invocation of method 'operator[]' ...
	// }

	return boost::report_errors();
}  // NOLINT(readability/fn_size)

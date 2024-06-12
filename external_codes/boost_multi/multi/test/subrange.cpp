// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <numeric>  // for std::iota

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

BOOST_AUTO_TEST_CASE(multi_array_range_section) {
	{
	#ifndef _MSC_VER
		multi::array<double, 4> arr({10, 20, 30, 40}, 99.0);
	#else
		multi::array<double, 4> arr(multi::extensions_t<4>{10, 20, 30, 40}, 99.0);
	#endif
		std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

		{
			static_assert(decltype(arr({0, 10}, {0, 20}, {0, 30}, {0, 40}))::rank::value == 4);
			static_assert(decltype(arr(5, {0, 20}, {0, 30}, {0, 40}))::rank::value == 3);
			static_assert(decltype(arr({0, 10}, 10, {0, 30}, {0, 40}))::rank::value == 3);
			static_assert(decltype(arr({0, 10}, {0, 20}, 15, {0, 40}))::rank::value == 3);
			static_assert(decltype(arr({0, 10}, {0, 20}, {0, 30}, 20))::rank::value == 3);

			static_assert(decltype(arr(5, 6, {0, 30}, {0, 40}))::rank::value == 2);
			static_assert(decltype(arr({0, 10}, 6, 15, {0, 40}))::rank::value == 2);
			static_assert(decltype(arr({0, 10}, {0, 20}, 15, 20))::rank::value == 2);

			static_assert(decltype(arr({0, 10}, {0, 20}, {0, 30}, {0, 40}))::rank_v == 4);
			static_assert(decltype(arr(5, {0, 20}, {0, 30}, {0, 40}))::rank_v == 3);
			static_assert(decltype(arr({0, 10}, 10, {0, 30}, {0, 40}))::rank_v == 3);
			static_assert(decltype(arr({0, 10}, {0, 20}, 15, {0, 40}))::rank_v == 3);
			static_assert(decltype(arr({0, 10}, {0, 20}, {0, 30}, 20))::rank_v == 3);

			static_assert(decltype(arr(5, 6, {0, 30}, {0, 40}))::rank_v == 2);
			static_assert(decltype(arr({0, 10}, 6, 15, {0, 40}))::rank_v == 2);
			static_assert(decltype(arr({0, 10}, {0, 20}, 15, 20))::rank_v == 2);
		}
		{
			auto&& all = arr({0, 10}, {0, 20}, {0, 30}, {0, 40});
			BOOST_REQUIRE( &arr[1][2][3][4] == &all[1][2][3][4] );
			BOOST_REQUIRE( &arr[1][2][3][4] == &arr({0, 10}, {0, 20}, {0, 30}, {0, 40})[1][2][3][4] );
		}
		{
			using multi::_;
			auto&& all = arr({0, 10}, {0, 20});
			BOOST_REQUIRE( &arr[1][2][3][4] == &all[1][2][3][4] );
		}
		{
			BOOST_REQUIRE( &arr(0, 0, 0, 0) == &arr[0][0][0][0] );
		}
		{
			auto&& sub = arr({0, 5}, {0, 10}, {0, 15}, {0, 20});
			BOOST_REQUIRE( &sub[1][2][3][4] == &arr[1][2][3][4] );
		}
	}
	{
		multi::array<double, 2> arr = {
			{1.0, 2.0, 3.0, 4.0},
			{5.0, 6.0, 7.0, 8.0},
			{9.0, 0.0, 1.0, 2.0},
			{3.0, 4.0, 5.0, 6.0},
		};
		multi::array<double, 2> arr2 = {
			{91.0, 92.0, 93.0, 94.0},
			{95.0, 96.0, 97.0, 98.0},
			{99.0, 90.0, 91.0, 92.0},
			{93.0, 94.0, 95.0, 96.0},
		};

		arr({0, 2}, {0, 2}) = arr2({0, 2}, {0, 2});
		BOOST_REQUIRE( arr != arr2 );
		BOOST_REQUIRE( arr({0, 2}, {0, 2}) == arr2({0, 2}, {0, 2}) );
		BOOST_REQUIRE( arr[1][1] == 96. );
	}
}

BOOST_AUTO_TEST_CASE(subrange_assignment) {
	multi::array<double, 2> const arr = {
		{1.0, 2.0, 3.0, 4.0},
		{5.0, 6.0, 7.0, 8.0},
		{9.0, 0.0, 1.0, 2.0},
		{3.0, 4.0, 5.0, 6.0},
	};
	{
		multi::array<double, 2> arr2 = {
			{9.0, 9.0, 9.0},
			{9.0, 9.0, 9.0},
			{9.0, 9.0, 9.0},
		};
		arr2({0, 3}, {0, 3}) = arr({0, 3}, {0, 3});
		BOOST_REQUIRE( arr2[1][2] == arr[1][2] );
	}
	{
		multi::array<double, 2> arr2 = {
			{9.0, 9.0, 9.0},
			{9.0, 9.0, 9.0},
			{9.0, 9.0, 9.0},
		};
		arr2() = arr({0, 3}, {0, 3});
		BOOST_REQUIRE( arr2[1][2] == arr[1][2] );
		BOOST_REQUIRE( arr2() == arr({0, 3}, {0, 3}) );
	}
	{
		multi::array<double, 2> arr2 = {
			{9.0, 9.0, 9.0},
			{9.0, 9.0, 9.0},
			{9.0, 9.0, 9.0},
		};
		arr2 = arr({0, 3}, {0, 3});
		BOOST_REQUIRE( arr2[1][2] == arr[1][2] );
		BOOST_REQUIRE( arr2 == arr({0, 3}, {0, 3}) );
	}
}

BOOST_AUTO_TEST_CASE(subrange_ranges_sliced_1D) {
	multi::array<double, 1> arr = {1.0, 2.0, 3.0, 4.0};
	auto&&                  Ab  = arr.sliced(1, 3);
	BOOST_REQUIRE( &Ab[0] == &arr[1] );

	auto&& Ab2 = Ab;
	BOOST_REQUIRE( &Ab2[0] == &arr[1] );

	//  auto Abb = Ab;  // not allowed!
	//  auto Abb = std::move(Ab); (void)Abb;

	auto const& Abc = arr.sliced(1, 3);
	BOOST_REQUIRE( &Abc[0] == &arr[1] );

	auto Aba = arr.sliced(1, 3);
	BOOST_REQUIRE( &Aba[0] == &arr[1] );
}

BOOST_AUTO_TEST_CASE(subrange_ranges_sliced) {
	multi::array<double, 2> arr = {
		{1.0, 2.0, 3.0, 4.0},
		{5.0, 6.0, 7.0, 8.0},
		{9.0, 0.0, 1.0, 2.0},
		{3.0, 4.0, 5.0, 6.0},
	};
	auto&& Ab = arr.sliced(0, 3);
	BOOST_REQUIRE( &Ab[2][2] == &arr[2][2] );

	auto const& Abc = arr.sliced(0, 3);
	BOOST_REQUIRE( &Abc[2][2] == &arr[2][2] );

	auto AB = arr.sliced(0, 3);
	BOOST_REQUIRE( &AB[2][2] == &arr[2][2] );
}

BOOST_AUTO_TEST_CASE(subrange_ranges) {
	multi::array<double, 2> arr = {
		{1.0, 2.0, 3.0, 4.0},
		{5.0, 6.0, 7.0, 8.0},
		{9.0, 0.0, 1.0, 2.0},
		{3.0, 4.0, 5.0, 6.0},
	};
	auto&& Ab = arr({0, 3}, {0, 3});
	BOOST_REQUIRE( &Ab[2][2] == &arr[2][2] );

	auto const& Abc = arr({0, 3}, {0, 3});
	BOOST_REQUIRE( &Abc[2][2] == &arr[2][2] );

	auto AB = arr({0, 3}, {0, 3});
	BOOST_REQUIRE( &AB[2][2] == &arr[2][2] );
}

BOOST_AUTO_TEST_CASE(subrange_1D_issue129) {
	multi::array<double, 1> arr({1024}, double{});
	std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

	BOOST_REQUIRE( arr.sliced(0, 512, 2)[  1] ==   2.0 );
	BOOST_REQUIRE( arr.sliced(0, 512, 2)[255] == 510.0 );

	BOOST_REQUIRE( arr.sliced(0, 512)[  1] ==   1.0 );
	BOOST_REQUIRE( arr.sliced(0, 512)[511] == 511.0 );

	BOOST_REQUIRE( arr({0, 512})[  1] ==   1.0 );
	BOOST_REQUIRE( arr({0, 512})[511] == 511.0 );

	//  BOOST_REQUIRE( arr({0, 512, 2})[  1] ==   2. );  // TODO(correaa) coompilation error
	//  BOOST_REQUIRE( arr({0, 512, 2})[255] == 510. );  // TODO(correaa) coompilation error
}

BOOST_AUTO_TEST_CASE(subrange_2D_issue129) {
	multi::array<double, 2> arr({66, 1024}, double{});
	std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

	BOOST_REQUIRE( arr[0].sliced(0, 512, 2)[  1] ==   2.0 );
	BOOST_REQUIRE( arr[0].sliced(0, 512, 2)[255] == 510.0 );

	BOOST_REQUIRE( arr[0].sliced(0, 512)[  1] ==   1.0 );
	BOOST_REQUIRE( arr[0].sliced(0, 512)[511] == 511.0 );

	BOOST_REQUIRE( arr(0, {0, 512})[  1] ==   1.0 );
	BOOST_REQUIRE( arr(0, {0, 512})[511] == 511.0 );

	//  BOOST_REQUIRE( arr(0, {0, 512, 2})[  1] ==   2. );  // TODO(correaa) coompilation error
	//  BOOST_REQUIRE( arr(0, {0, 512, 2})[255] == 510. );  // TODO(correaa) coompilation error
}

class rng3_t {
	int start_;
	int finish_;

 public:
	rng3_t(int start, int finish) : start_{start}, finish_{finish} {}  // NOLINT(bugprone-easily-swappable-parameters)
	auto first() const { return start_; }
	auto last() const { return finish_; }
};

BOOST_AUTO_TEST_CASE(subrange_start_finish) {
	multi::array<double, 2> arr = {
		{ 1.0,  2.0},
		{ 3.0,  4.0},
		{ 5.0,  6.0},
		{ 7.0,  8.0},
		{ 9.0, 10.0},
		{11.0, 12.0},
		{13.0, 14.0},
	};
	BOOST_REQUIRE( &arr({2, 5}, 1)[0] == &arr[2][1] );

	multi::irange const rng(2, 5);
	BOOST_REQUIRE( &arr(rng, 1)[0] == &arr[2][1] );

	struct : multi::irange {
		using multi::irange::irange;
	} const rng2(2, 5);

	BOOST_REQUIRE( &arr(rng2, 1)[0] == &arr[2][1] );

	rng3_t const rng3{2, 5};

	multi::irange const rng4(rng3);

	BOOST_REQUIRE( &arr(rng4, 1)[0] == &arr[2][1] );

	BOOST_REQUIRE( &arr(rng3, 1)[0] == &arr[2][1] );
}

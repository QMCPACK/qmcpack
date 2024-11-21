// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <numeric>  // for std::iota

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE)  /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(multi_array_range_section) {
	{
#ifndef _MSC_VER
		multi::array<double, 4> arr({ 10, 20, 30, 40 }, 99.0);
#else
		multi::array<double, 4> arr(multi::extensions_t<4>{ 10, 20, 30, 40 }, 99.0);
#endif
		std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

		{
			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, { 0, 30 }, { 0, 40 }))::rank::value == 4);
			static_assert(decltype(arr(5, { 0, 20 }, { 0, 30 }, { 0, 40 }))::rank::value == 3);
			static_assert(decltype(arr({ 0, 10 }, 10, { 0, 30 }, { 0, 40 }))::rank::value == 3);
			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, 15, { 0, 40 }))::rank::value == 3);
			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, { 0, 30 }, 20))::rank::value == 3);

			static_assert(decltype(arr(5, 6, { 0, 30 }, { 0, 40 }))::rank::value == 2);
			static_assert(decltype(arr({ 0, 10 }, 6, 15, { 0, 40 }))::rank::value == 2);
			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, 15, 20))::rank::value == 2);

			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, { 0, 30 }, { 0, 40 }))::rank_v == 4);
			static_assert(decltype(arr(5, { 0, 20 }, { 0, 30 }, { 0, 40 }))::rank_v == 3);
			static_assert(decltype(arr({ 0, 10 }, 10, { 0, 30 }, { 0, 40 }))::rank_v == 3);
			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, 15, { 0, 40 }))::rank_v == 3);
			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, { 0, 30 }, 20))::rank_v == 3);

			static_assert(decltype(arr(5, 6, { 0, 30 }, { 0, 40 }))::rank_v == 2);
			static_assert(decltype(arr({ 0, 10 }, 6, 15, { 0, 40 }))::rank_v == 2);
			static_assert(decltype(arr({ 0, 10 }, { 0, 20 }, 15, 20))::rank_v == 2);
		}
		{
			auto&& all = arr({ 0, 10 }, { 0, 20 }, { 0, 30 }, { 0, 40 });
			BOOST_TEST( &arr[1][2][3][4] == &all[1][2][3][4] );
			BOOST_TEST( &arr[1][2][3][4] == &arr({0, 10}, {0, 20}, {0, 30}, {0, 40})[1][2][3][4] );
		}
		{
			using multi::_;
			auto&& all = arr({ 0, 10 }, { 0, 20 });
			BOOST_TEST( &arr[1][2][3][4] == &all[1][2][3][4] );
		}
		{
			BOOST_TEST( &arr(0, 0, 0, 0) == &arr[0][0][0][0] );
		}
		{
			auto&& sub = arr({ 0, 5 }, { 0, 10 }, { 0, 15 }, { 0, 20 });
			BOOST_TEST( &sub[1][2][3][4] == &arr[1][2][3][4] );
		}
	}
	{
		multi::array<int, 2> arr = {
			{ 10, 20, 30, 40 },
			{ 50, 60, 70, 80 },
			{ 90, 00, 10, 20 },
			{ 30, 40, 50, 60 },
		};
		multi::array<int, 2> arr2 = {
			{ 910, 920, 930, 940 },
			{ 950, 960, 970, 980 },
			{ 990, 900, 910, 920 },
			{ 930, 940, 950, 960 },
		};

		arr({ 0, 2 }, { 0, 2 }) = arr2({ 0, 2 }, { 0, 2 });
		BOOST_TEST( arr != arr2 );
		BOOST_TEST( arr({0, 2}, {0, 2}) == arr2({0, 2}, {0, 2}) );
		BOOST_TEST( arr[1][1] == 960 );
	}
}

BOOST_AUTO_TEST_CASE(subrange_assignment) {
	multi::array<int, 2> const arr = {
		{ 10, 20, 30, 40 },
		{ 50, 60, 70, 80 },
		{ 90, 00, 10, 20 },
		{ 30, 40, 50, 60 },
	};
	{
		multi::array<int, 2> arr2 = {
			{ 90, 90, 90 },
			{ 90, 90, 90 },
			{ 90, 90, 90 },
		};
		arr2({ 0, 3 }, { 0, 3 }) = arr({ 0, 3 }, { 0, 3 });
		BOOST_TEST( arr2[1][2] == arr[1][2] );
	}
	{
		multi::array<int, 2> arr2 = {
			{ 90, 90, 90 },
			{ 90, 90, 90 },
			{ 90, 90, 90 },
		};
		arr2() = arr({ 0, 3 }, { 0, 3 });
		BOOST_TEST( arr2[1][2] == arr[1][2] );
		BOOST_TEST( arr2() == arr({0, 3}, {0, 3}) );
	}
	{
		multi::array<int, 2> arr2 = {
			{ 90, 90, 90 },
			{ 90, 90, 90 },
			{ 90, 90, 90 },
		};
		arr2 = arr({ 0, 3 }, { 0, 3 });
		BOOST_TEST( arr2[1][2] == arr[1][2] );
		BOOST_TEST( arr2 == arr({0, 3}, {0, 3}) );
	}
}

BOOST_AUTO_TEST_CASE(subrange_ranges_sliced_1D) {
	multi::array<double, 1> arr = { 1.0, 2.0, 3.0, 4.0 };
	auto&&                  Ab  = arr.sliced(1, 3);
	BOOST_TEST( &Ab[0] == &arr[1] );

	auto&& Ab2 = Ab;
	BOOST_TEST( &Ab2[0] == &arr[1] );

	//  auto Abb = Ab;  // not allowed!
	//  auto Abb = std::move(Ab); (void)Abb;

	auto const& Abc = arr.sliced(1, 3);
	BOOST_TEST( &Abc[0] == &arr[1] );

	auto Aba = arr.sliced(1, 3);
	BOOST_TEST( &Aba[0] == &arr[1] );
}

BOOST_AUTO_TEST_CASE(subrange_ranges_sliced) {
	multi::array<double, 2> arr = {
		{ 1.0, 2.0, 3.0, 4.0 },
		{ 5.0, 6.0, 7.0, 8.0 },
		{ 9.0, 0.0, 1.0, 2.0 },
		{ 3.0, 4.0, 5.0, 6.0 },
	};
	auto&& Ab = arr.sliced(0, 3);
	BOOST_TEST( &Ab[2][2] == &arr[2][2] );

	auto const& Abc = arr.sliced(0, 3);
	BOOST_TEST( &Abc[2][2] == &arr[2][2] );

	auto AB = arr.sliced(0, 3);
	BOOST_TEST( &AB[2][2] == &arr[2][2] );
}

BOOST_AUTO_TEST_CASE(subrange_ranges) {
	multi::array<double, 2> arr = {
		{ 1.0, 2.0, 3.0, 4.0 },
		{ 5.0, 6.0, 7.0, 8.0 },
		{ 9.0, 0.0, 1.0, 2.0 },
		{ 3.0, 4.0, 5.0, 6.0 },
	};
	auto&& Ab = arr({ 0, 3 }, { 0, 3 });
	BOOST_TEST( &Ab[2][2] == &arr[2][2] );

	auto const& Abc = arr({ 0, 3 }, { 0, 3 });
	BOOST_TEST( &Abc[2][2] == &arr[2][2] );

	auto AB = arr({ 0, 3 }, { 0, 3 });
	BOOST_TEST( &AB[2][2] == &arr[2][2] );
}

BOOST_AUTO_TEST_CASE(subrange_1D_issue129) {
	multi::array<int, 1> arr({ 1024 }, int{});
	std::iota(arr.elements().begin(), arr.elements().end(), 0);

	BOOST_TEST( arr.sliced(0, 512, 2)[  1] ==   2 );
	BOOST_TEST( arr.sliced(0, 512, 2)[255] == 510 );

	BOOST_TEST( arr.sliced(0, 512)[  1] ==   1 );
	BOOST_TEST( arr.sliced(0, 512)[511] == 511 );

	BOOST_TEST( arr({0, 512})[  1] ==   1 );
	BOOST_TEST( arr({0, 512})[511] == 511 );

	//  BOOST_TEST( arr({0, 512, 2})[  1] ==   2 );  // TODO(correaa) coompilation error
	//  BOOST_TEST( arr({0, 512, 2})[255] == 510 );  // TODO(correaa) coompilation error
}

BOOST_AUTO_TEST_CASE(subrange_2D_issue129) {
	multi::array<int, 2> arr({ 66, 1024 }, int{});
	std::iota(arr.elements().begin(), arr.elements().end(), 0);

	BOOST_TEST( arr[0].sliced(0, 512, 2)[  1] ==   2 );
	BOOST_TEST( arr[0].sliced(0, 512, 2)[255] == 510 );

	BOOST_TEST( arr[0].sliced(0, 512)[  1] ==   1 );
	BOOST_TEST( arr[0].sliced(0, 512)[511] == 511 );

	BOOST_TEST( arr(0, {0, 512})[  1] ==   1 );
	BOOST_TEST( arr(0, {0, 512})[511] == 511 );

	// BOOST_TEST( arr(0, {0, 512, 2})[  1] ==   2 );  // TODO(correaa) coompilation error
	// BOOST_TEST( arr(0, {0, 512, 2})[255] == 510 );  // TODO(correaa) coompilation error
}

class rng3_t {
	int start_;
	int finish_;

 public:
	rng3_t(int start, int finish) : start_{ start }, finish_{ finish } {}  // NOLINT(bugprone-easily-swappable-parameters)
	auto first() const { return start_; }
	auto last() const { return finish_; }
};

BOOST_AUTO_TEST_CASE(subrange_start_finish) {
	multi::array<double, 2> arr = {
		{  1.0,  2.0 },
		{  3.0,  4.0 },
		{  5.0,  6.0 },
		{  7.0,  8.0 },
		{  9.0, 10.0 },
		{ 11.0, 12.0 },
		{ 13.0, 14.0 },
	};
	BOOST_TEST( &arr({2, 5}, 1)[0] == &arr[2][1] );

	multi::irange const rng(2, 5);
	BOOST_TEST( &arr(rng, 1)[0] == &arr[2][1] );

	struct : multi::irange {
		using multi::irange::irange;
	} const rng2(2, 5);

	BOOST_TEST( &arr(rng2, 1)[0] == &arr[2][1] );

	rng3_t const rng3{ 2, 5 };

	multi::irange const rng4(rng3);

	BOOST_TEST( &arr(rng4, 1)[0] == &arr[2][1] );

	BOOST_TEST( &arr(rng3, 1)[0] == &arr[2][1] );
}
return boost::report_errors();}

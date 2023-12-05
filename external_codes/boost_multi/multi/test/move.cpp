// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi move"
#include<boost/test/unit_test.hpp>

#include <boost/multi_array.hpp>

#include "multi/array.hpp"

#include <algorithm>  // for std::move
#include <memory>
#include <vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(move_unique_ptr_1D) {
	{
		multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
		arr[1] = std::make_unique<int>(42);

		multi::array<std::unique_ptr<int>, 1> arr2(multi::extensions_t<1>{10});
		std::move(arr.begin(), arr.end(), arr2.begin());

		BOOST_REQUIRE( !arr[1] );
		BOOST_REQUIRE(  arr2[1] );
		BOOST_REQUIRE( *arr2[1] == 42 );
	}
	{
		multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
		arr[1] = std::make_unique<int>(42);

		multi::array<std::unique_ptr<int>, 1> arr2 = std::move(arr);
		BOOST_REQUIRE(  arr.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)
		BOOST_REQUIRE(  arr2[1] );
		BOOST_REQUIRE( *arr2[1] == 42 );
	}
	{
		multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
		arr[1] = std::make_unique<int>(42);

		multi::array<std::unique_ptr<int>, 1> arr2;//(multi::extensions_t<1>{10});
		arr2 = std::move(arr);
		BOOST_REQUIRE(  arr.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)
		BOOST_REQUIRE(  arr2[1] );
		BOOST_REQUIRE( *arr2[1] == 42 );
	}
	{
		multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
		arr[1] = std::make_unique<int>(42);

		multi::array<std::unique_ptr<int>, 1> arr2(multi::extensions_t<1>{10});
	//  arr2() = arr();  // fails to compile, elements are not copy assignable
		arr2() = arr().moved();
		BOOST_REQUIRE( !arr[1] );
		BOOST_REQUIRE(  arr2[1] );
		BOOST_REQUIRE( *arr2[1] == 42 );
	}
}

BOOST_AUTO_TEST_CASE(multi_swap) {
	multi::array<double, 2> arr({3,  5}, 99.);
	multi::array<double, 2> arr2({7, 11}, 88.);
	swap(arr, arr2);
	BOOST_REQUIRE( size(arr) == 7 );
	BOOST_REQUIRE( arr[1][2] == 88. );
	BOOST_REQUIRE( arr2[1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_std_swap) {
	multi::array<double, 2> arr({3,  5}, 99.);
	multi::array<double, 2> arr2({7, 11}, 88.);
	using std::swap;
	swap(arr, arr2);
	BOOST_REQUIRE( size(arr) == 7 );
	BOOST_REQUIRE( arr[1][2] == 88. );
	BOOST_REQUIRE( arr2[1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_clear) {
	multi::array<double, 2> arr({10, 10}, 99.);
	arr.clear();
	BOOST_REQUIRE(arr.is_empty());
	arr.reextent({20, 20}, 99.);
	BOOST_REQUIRE(not arr.is_empty());
	clear(arr).reextent({30, 30}, 88.);
	BOOST_REQUIRE(arr[15][15] == 88.);
}

BOOST_AUTO_TEST_CASE(multi_array_move) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	multi::array<double, 2> arr2(std::move(Av[0]), std::allocator<double>{});

	BOOST_REQUIRE( is_empty(Av[0]) );
	BOOST_REQUIRE( size(arr2) == 4 );
	BOOST_REQUIRE( arr2[1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv; Bv.reserve(Av.size());

	std::move( begin(Av), end(Av), std::back_inserter(Bv) );

	BOOST_REQUIRE( size(Bv) == size(Av) );
	BOOST_REQUIRE( is_empty(Av[4]) );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_reserve) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv; Bv.reserve(Av.size());

//	for(auto& v: Av) Bv.emplace_back(std::move(v), std::allocator<double>{}); // segfaults nvcc 11.0 but not nvcc 11.1
	std::move(begin(Av), end(Av), std::back_inserter(Bv));

	BOOST_REQUIRE( size(Bv) == size(Av) );
	BOOST_REQUIRE( is_empty(Av[4]) );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_move) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv = std::move(Av);

	Av.clear();
	BOOST_REQUIRE( size(Av) == 0 );
	BOOST_REQUIRE( size(Bv) == 10 );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_array) {
	multi::array<std::vector<double>, 2> arr({10, 10}, std::vector<double>(5) );
	auto arr2 = std::move(arr);
	BOOST_REQUIRE( arr .   empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move) test deterministic moved from state
	BOOST_REQUIRE( arr .is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move) test deterministic moved from state
	BOOST_REQUIRE( arr2.size() == 10 );
}

BOOST_AUTO_TEST_CASE(multi_array_move_elements) {
	multi::array<std::vector<double>, 1> arr({10}, std::vector<double>(5) );

	std::vector<std::vector<double>> sink(5);

	auto* ptr1 = arr[1].data();

	std::copy( arr({0, 5}).moved().begin(), arr({0, 5}).moved().end(), sink.begin() );
	BOOST_REQUIRE(     arr[1].empty() );
	BOOST_REQUIRE( not arr[5].empty() );

	BOOST_REQUIRE( sink[1].data() == ptr1 );
}

BOOST_AUTO_TEST_CASE(multi_array_move_elements_range) {
	multi::array<std::vector<double>, 1> arr({10}, std::vector<double>(5) );

	std::vector<std::vector<double>> sink(5);

	auto* ptr1 = arr[1].data();

	std::copy( arr({0, 5}).moved().elements().begin(), arr({0, 5}).moved().elements().end(), sink.begin() );
	BOOST_REQUIRE(     arr[1].empty() );
	BOOST_REQUIRE( not arr[5].empty() );

	BOOST_REQUIRE( sink[1].data() == ptr1 );
}

BOOST_AUTO_TEST_CASE(multi_array_move_elements_to_array) {
	multi::array<std::vector<double>, 1> arr({10}, std::vector<double>(5, 99.) );
	BOOST_REQUIRE( arr.size() == 10 );
	multi::array<std::vector<double>, 1> arr2({ 5}, {}, {});

	auto* ptr1 = arr[1].data();

	arr2().elements() = arr({0, 5}).moved().elements();

	BOOST_REQUIRE( arr2[1].size() == 5 );
	BOOST_REQUIRE( arr2[1][4] == 99. );

	BOOST_REQUIRE(     arr[1].empty() );
	BOOST_REQUIRE( not arr[5].empty() );

	BOOST_REQUIRE( arr2[1].data() == ptr1 );
}

BOOST_AUTO_TEST_CASE(move_range_vector_1D) {
	std::vector<std::vector<double>> arr(10, std::vector<double>{1., 2., 3.});
	std::vector<std::vector<double>> arr2(10);
	std::move(arr.begin(), arr.end(), arr2.begin());

	BOOST_REQUIRE( arr2[0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[0].empty() );
	BOOST_REQUIRE( arr[1].empty() );
}

BOOST_AUTO_TEST_CASE(copy_range_1D) {
	multi::array<std::vector<double>, 1> arr({3}, std::vector<double>{1., 2., 3.});
	BOOST_REQUIRE( arr.size() == 3 );
	multi::array<std::vector<double>, 1> arr2({3}, std::vector<double>{});
	std::copy(arr.begin(), arr.end(), arr2.begin());

	BOOST_REQUIRE( arr2[0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr[1] == std::vector<double>({1., 2., 3.}) );
}

BOOST_AUTO_TEST_CASE(move_range_1D) {
	multi::array<std::vector<double>, 1> arr({3}, std::vector<double>{1., 2., 3.});
	BOOST_REQUIRE( arr.size() == 3 );
	multi::array<std::vector<double>, 1> arr2({3}, std::vector<double>{});
	std::move(arr.begin(), arr.end(), arr2.begin());

	BOOST_REQUIRE( arr2[0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[0].empty() );
	BOOST_REQUIRE( arr[1].empty() );
}

BOOST_AUTO_TEST_CASE(move_range_1D_moved_begin) {
	multi::array<std::vector<double>, 1> arr({3}, std::vector<double>{1., 2., 3.});
	BOOST_REQUIRE( arr.size() == 3 );
	multi::array<std::vector<double>, 1> arr2({3}, std::vector<double>{});
	std::copy(arr.mbegin(), arr.mend(), arr2.begin());

	BOOST_REQUIRE( arr2[0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[0].empty() );
	BOOST_REQUIRE( arr[1].empty() );
}

template<class... Ts> void what(Ts&&...) = delete;

BOOST_AUTO_TEST_CASE(copy_move_range) {
	multi::array<std::vector<double>, 2> arr({10, 20}, std::vector<double>{1., 2., 3.});
	multi::array<std::vector<double>, 2> arr2({10, 20}, std::vector<double>{}          );

	std::copy(arr.mbegin(), arr.mend(), arr2.begin());

	BOOST_REQUIRE( arr2[0][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[0][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr2[1][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[0][0].empty() );
	BOOST_REQUIRE( arr[0][1].empty() );

	BOOST_REQUIRE( arr[1][0].empty() );
	BOOST_REQUIRE( arr[1][1].empty() );
}

BOOST_AUTO_TEST_CASE(copy_move_range_moved_begin) {
	multi::array<std::vector<double>, 2> arr({10, 20}, std::vector<double>{1., 2., 3.});
	multi::array<std::vector<double>, 2> arr2({10, 20}, std::vector<double>{}          );

	std::copy(arr.moved().begin(), arr.moved().end(), arr2.begin());

	BOOST_REQUIRE( arr2[0][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[0][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr2[1][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[0][0].empty() );
	BOOST_REQUIRE( arr[0][1].empty() );

	BOOST_REQUIRE( arr[1][0].empty() );
	BOOST_REQUIRE( arr[1][1].empty() );
}

BOOST_AUTO_TEST_CASE(copy_move_range_moved_begin_block) {
	multi::array<std::vector<double>, 2> arr({10, 20}, std::vector<double>{1., 2., 3.});
	multi::array<std::vector<double>, 2> arr2({ 3,  5}, std::vector<double>{}          );

	std::copy(arr({5, 8}, {10, 15}).moved().begin(), arr({5, 8}, {10, 15}).moved().end(), arr2.begin());

	BOOST_REQUIRE( arr2[0][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[0][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr2[1][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[5][10].empty() );
	BOOST_REQUIRE( arr[5][11].empty() );

	BOOST_REQUIRE( arr[6][10].empty() );
	BOOST_REQUIRE( arr[6][11].empty() );
}


BOOST_AUTO_TEST_CASE(move_reference_range) {
	multi::array<std::vector<double>, 2> arr({10, 20}, std::vector<double>{1., 2., 3.});
	multi::array<std::vector<double>, 2> arr2({10, 20}, std::vector<double>{}          );

//  arr2() = arr().moved();
	std::copy(arr().moved().begin(), arr().moved().end(), arr2().begin());

	BOOST_REQUIRE( arr2[0][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[0][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr2[1][0] == std::vector<double>({1., 2., 3.}) );
	BOOST_REQUIRE( arr2[1][1] == std::vector<double>({1., 2., 3.}) );

	BOOST_REQUIRE( arr[0][0].empty() );
	BOOST_REQUIRE( arr[0][1].empty() );

	BOOST_REQUIRE( arr[1][0].empty() );
	BOOST_REQUIRE( arr[1][1].empty() );
}

BOOST_AUTO_TEST_CASE(move_array_elements) {  // NOLINT(readability-function-cognitive-complexity)
	{
		auto arr = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = std::move(arr);
		BOOST_REQUIRE( arr2.size() == 5 );
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)
	}
	{
		auto arr = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));

		std::vector<double> v0 = std::move(arr[0]);
		BOOST_REQUIRE( v0.size() == 7 );
		BOOST_REQUIRE( arr[0].empty() );

		std::vector<double> v1 = std::move(arr)[1];
		BOOST_REQUIRE( v1.size() == 7 );
		BOOST_REQUIRE( arr[1].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move) for test

		auto arr2 = multi::array<std::vector<double>, 1>({ 1}, std::vector<double>{});

		arr2({0, 1}) = arr({2, 3});
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[2].size() == 7 );
	}
	{
		auto arr = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});

		arr2() = arr();
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[2].size() == 7 );

		arr2() = std::move(arr)();
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[2].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});

		arr2() = arr();
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].size() == 7 );
	}
	{
		auto arr = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});

		arr2() = std::move(arr)();
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});

		auto&& mAp = std::move(arr)();
		arr2() = mAp;
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});

		arr2({0, 5}) = std::move(arr)();
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});

		arr2() = arr.take(5);
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].size() == 7);  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});

		arr2() = std::move(arr).take(5);
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});
		auto&& mAt5 = std::move(arr).take(5);
		arr2() = mAt5;
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});
		auto&& mAt5 = std::move(arr).take(5);
		arr2() = mAt5;
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});
		auto&& mAt5 = std::move(arr).take(5);
		arr2() = std::move(mAt5);
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});
		auto&& mAt5 = std::move(arr).take(5);
		arr2() = std::move(mAt5).take(5);
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});
		auto&& mAt5 = std::move(arr).take(5);
		auto&& mAt5t5 = std::move(mAt5).take(5);
		arr2() = mAt5t5;
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>{});
		auto&& mAt5 = std::move(arr).take(5);
		arr2() = std::move(mAt5).drop(0);
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
	{
		auto arr  = multi::array<std::vector<double>, 1>({ 5}, std::vector<double>(7));
		auto arr2 = multi::array<std::vector<double>, 1>({ 4}, std::vector<double>{});
		arr2() = std::move(arr).drop(1);
		BOOST_REQUIRE( arr2[0].size() == 7 );
		BOOST_REQUIRE( arr[0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		BOOST_REQUIRE( arr[1].empty()     );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	}
}

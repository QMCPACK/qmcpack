// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <array>
#include <deque>
#include <numeric>  // for std::iota

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
// #  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
// #  pragma GCC diagnostic ignored "-Wfloat-equal"
#elif defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4244)
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(empty_intersection) {
	multi::array<double, 1> arr({10});
	multi::array<double, 1> arr2;

	auto const is = intersection(arr.extension(), arr2.extension());
	BOOST_REQUIRE( arr(is).is_empty() );
	arr2(is) = arr(is);

	BOOST_REQUIRE( arr2(is) == arr(is) );
}

BOOST_AUTO_TEST_CASE(multi_tests_element_access_with_tuple) {
	multi::array<char, 2> arr({3, 3}, 'k');

	std::array<int, 2> point = {
		{1, 2}
	};

	BOOST_REQUIRE(  arr[point[0]][point[1]] ==  arr(1, 2) );
	BOOST_REQUIRE( &arr(point[0], point[1]) == &arr[point[0]][point[1]] );

	BOOST_REQUIRE( &arr[point[0]][point[1]] == &arr(point[0], point[1]) );
	BOOST_REQUIRE( &arr(point[0], point[1]) == &arr.apply(point) );

	BOOST_REQUIRE( &arr[point[0]][point[1]] == &std::apply(arr, point) );
	BOOST_REQUIRE( &arr[point[0]][point[1]] == &     apply(arr, point) );
}

BOOST_AUTO_TEST_CASE(multi_tests_extension_with_tuple) {
	{
		multi::array<double, 2>::extensions_type const ext = {3, 4};

		multi::array<double, 2> const arr(ext, 44.0);

		BOOST_REQUIRE( size(arr) == 3 );
	}
	{
		auto const [en, em] = std::make_tuple(3, 4);
		multi::array<double, 2> const arr({en, em}, 44.0);
		BOOST_REQUIRE( size(arr) == 3 );
	}
	{
		auto arr = std::apply([](auto const&... szs) { return multi::array<double, 2>({szs...}, 55.0); }, std::make_tuple(3, 4));
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( std::get<0>(sizes(arr)) == 3 );
		BOOST_REQUIRE( std::get<1>(sizes(arr)) == 4 );
	}
}

BOOST_AUTO_TEST_CASE(multi_test_constness_reference) {
	multi::array<char, 2> const carr({10, 10}, '9');

	BOOST_REQUIRE( size( carr(1, {0, 3}) ) == 3 );

	BOOST_REQUIRE( carr(1, {0, 3})[1] == '9' );
	static_assert(decltype(carr({0, 3}, 1))::rank_v == 1);
	BOOST_REQUIRE( size(carr.sliced(0, 3)) == 3 );

	BOOST_REQUIRE( carr.range({0, 3}).rotated()[1].unrotated().size() == 3 );

	BOOST_REQUIRE( carr({0, 3}, {0, 3})[1][1] == '9' );

	static_assert(! std::is_assignable_v<decltype(carr(1, {0, 3})[1]), double>);
}

BOOST_AUTO_TEST_CASE(multi_test_stencil) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) ""s

	multi::array<std::string, 2> arr = {
		{"a"s, "b"s, "c"s, "d"s, "e"s},
		{"f"s, "g"s, "h"s, "f"s, "g"s},
		{"h"s, "i"s, "j"s, "k"s, "l"s},
	};

	BOOST_REQUIRE(      size(arr) == 3                                            );
	BOOST_REQUIRE(           arr.num_elements() == 3*5L                           );
	BOOST_REQUIRE(           arr[1][2] == "h"                                     );

	BOOST_REQUIRE(      size(arr          ({1, 3}, {2, 5})) == 2                  );
	BOOST_REQUIRE( extension(arr          ({1, 3}, {2, 5})).first() == 0          );
	BOOST_REQUIRE(           arr          ({1, 3}, {2, 5}).num_elements() == 2*3L );
	BOOST_REQUIRE(           arr          ({1, 3}, {2, 5}).num_elements() == 2*3L );
	BOOST_REQUIRE(           arr          ({1, 3}, {2, 5})[0][0] == "h"           );
	BOOST_REQUIRE(          &arr          ({1, 3}, {2, 5})[0][0] == &arr[1][2]      );

	BOOST_REQUIRE(      size(arr.stenciled({1, 3}, {2, 5})) == 2                  );
	BOOST_REQUIRE( extension(arr.stenciled({1, 3}, {2, 5})).first() == 1          );
	BOOST_REQUIRE(           arr.stenciled({1, 3}, {2, 5}).num_elements() == 2*3L );
	BOOST_REQUIRE(           arr.stenciled({1, 3}, {2, 5}) [1][2] == "h"          );
	BOOST_REQUIRE(          &arr.stenciled({1, 3}, {2, 5}) [1][2] == &arr[1][2]     );

	BOOST_REQUIRE(  arr().elements().size() == arr.num_elements() );

	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements()[0] == &arr(1, 2) );
	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements()[arr({1, 3}, {2, 5}).elements().size() - 1] == &arr(2, 4) );

	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements().front() == &arr(1, 2) );
	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements().back()  == &arr(2, 4) );
}

BOOST_AUTO_TEST_CASE(empty_elements) {
	multi::array<int, 2> arr1;
	multi::array<int, 2> arr2;

	BOOST_REQUIRE( arr1.elements().size() == 0 );
	BOOST_REQUIRE( arr2.elements().size() == 0 );
	BOOST_REQUIRE(   arr1.elements() == arr2.elements()  );
	BOOST_REQUIRE( !(arr1.elements() != arr2.elements()) );
}

BOOST_AUTO_TEST_CASE(multi_test_elements_1D) {
	multi::array<int, 1> arr = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	BOOST_REQUIRE( arr.size() == 10 );

	BOOST_REQUIRE(  arr.elements().size() == 10 );
	BOOST_REQUIRE( &arr.elements()[0] == &arr[0] );
	BOOST_REQUIRE( &arr.elements()[9] == &arr[9] );

	BOOST_REQUIRE(    arr.elements().begin() <  arr.elements().end()     );
	BOOST_REQUIRE(    arr.elements().end()   >  arr.elements().begin()   );
	BOOST_REQUIRE(    arr.elements().begin() != arr.elements().end()     );
	BOOST_REQUIRE( !( arr.elements().begin() == arr.elements().end()   ) );

	BOOST_REQUIRE(  arr().elements().begin() <  arr().elements().end() );
	BOOST_REQUIRE(  arr().elements().begin() == arr().elements().begin() );

	BOOST_REQUIRE( arr().elements().begin() <  arr().elements().end() || arr().elements().begin() == arr().elements().end() );
	BOOST_REQUIRE( arr().elements().begin() <= arr().elements().end() );

	BOOST_REQUIRE(  arr().elements().end()  >  arr().elements().begin() );
	BOOST_REQUIRE(  arr().elements().end()  >= arr().elements().begin() );

	arr.elements() = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
	BOOST_REQUIRE( arr[2] == 7 );
	BOOST_REQUIRE( arr.elements()[2] == 7 );
}

BOOST_AUTO_TEST_CASE(multi_test_elements_1D_as_range) {
	multi::array<int, 1> arr = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	BOOST_REQUIRE( arr.size() == 10 );

	arr().elements() = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
	BOOST_REQUIRE( arr[2] == 7 );
	BOOST_REQUIRE( arr.elements()[2] == 7 );
}

BOOST_AUTO_TEST_CASE(elements_from_init_list_2D) {
	multi::array<int, 2> arr({3, 2});
	arr().elements() = {1, 2, 3, 4, 5, 6};
	BOOST_REQUIRE(arr[1][0] == 3);

	arr.elements() = {10, 20, 30, 40, 50, 60};
	BOOST_REQUIRE(arr[1][0] == 30);
}

BOOST_AUTO_TEST_CASE(front_back_2D) {
	multi::array<int, 2> arr({3, 4});
	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 0);

	BOOST_REQUIRE(  arr.front()[2] ==  arr[0][2] );
	BOOST_REQUIRE( &arr.front()[2] == &arr[0][2] );

	BOOST_REQUIRE(  arr.back ()[2] ==  arr[2][2] );
	BOOST_REQUIRE( &arr.back ()[2] == &arr[2][2] );
}

BOOST_AUTO_TEST_CASE(front_back_1D) {
	multi::array<int, 1> arr({30}, double{});
	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 0);

	BOOST_REQUIRE(  arr.front() ==  arr[ 0] );
	BOOST_REQUIRE( &arr.front() == &arr[ 0] );

	BOOST_REQUIRE(  arr.back () ==  arr[29] );
	BOOST_REQUIRE( &arr.back () == &arr[29] );
}

BOOST_AUTO_TEST_CASE(elements_rvalues) {
	using movable_type = std::vector<int>;
	movable_type const movable_value(5, 99);  // NOLINT(fuchsia-default-arguments-calls)

	multi::array<movable_type, 1> arr = {movable_value, movable_value, movable_value};
	BOOST_REQUIRE( arr.size() == 3 );

	movable_type const front = std::move(arr)[0];

	BOOST_REQUIRE( front == movable_value );
	BOOST_REQUIRE( arr[0].empty()           );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes
	BOOST_REQUIRE( arr[1] == movable_value  );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

	std::move(arr)[1] = movable_value;
}

template<class Array1D>
void assign_elements_from_to(Array1D&& arr, std::deque<std::vector<double>>& dest) {  // NOLINT(google-runtime-references) dest is mutated
	std::copy(std::forward<Array1D>(arr).begin(), std::forward<Array1D>(arr).end(), std::back_inserter(dest));
}

BOOST_AUTO_TEST_CASE(elements_rvalues_nomove) {
	using movable_type = std::vector<double>;
	movable_type const movable_value(5., 99.0);  // NOLINT(fuchsia-default-arguments-calls)

	multi::array<movable_type, 1> arr = {movable_value, movable_value, movable_value};
	BOOST_REQUIRE( arr.size() == 3 );

	std::deque<std::vector<double>> q1;

	assign_elements_from_to(arr, q1);

	BOOST_REQUIRE( arr[0] == movable_value );

	std::deque<std::vector<double>> q2;

	assign_elements_from_to(std::move(arr), q2);

	//  BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

	BOOST_REQUIRE( q1 == q2 );
}

BOOST_AUTO_TEST_CASE(elements_rvalues_assignment) {
	std::vector<double> vec = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)

	std::move(vec) = std::vector<double>{3.0, 4.0, 5.0};  // NOLINT(fuchsia-default-arguments-calls)

	std::move(vec)[1] = 99.0;  // it compiles  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

	multi::array<double, 1>       arr1 = {1.0, 2.0, 3.0};
	multi::array<double, 1> const arr2 = {1.0, 2.0, 3.0};

	std::move(arr1) = arr2;  // this compiles TODO(correaa) should it?
}

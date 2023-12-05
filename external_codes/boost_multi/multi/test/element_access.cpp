// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi element access"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include <deque>
#include <numeric>  // for iota

namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

namespace test_bee {
	struct bee{};

	template<class Array> auto paren(Array&& arr, bee const&/*unused*/) -> decltype(auto) {
		return std::forward<Array>(arr)(0);
	}
}  // end namespace test_bee

BOOST_AUTO_TEST_CASE(overload_paren) {
	multi::array<double, 1> arr({10});
	test_bee::bee zero;
	BOOST_REQUIRE( &arr(0) == &arr(zero) );
}

BOOST_AUTO_TEST_CASE(empty_intersection) {
	multi::array<double, 1> arr({10});
	multi::array<double, 1> arr2;

	auto const is = intersection(arr.extension(), arr2.extension());
	BOOST_REQUIRE( arr(is).is_empty() );
	arr2(is) = arr(is);

	BOOST_REQUIRE( arr2(is) == arr(is) );
}

BOOST_AUTO_TEST_CASE(multi_tests_element_access_with_tuple) {
	multi::array<double, 2> arr({3, 3}, 44.);
	std::array<int, 2> point = {{1, 2}};

	BOOST_REQUIRE(  arr[point[0]][point[1]] ==  arr(1, 2) );
	BOOST_REQUIRE( &arr(point[0], point[1]) == &arr[point[0]][point[1]] );

	BOOST_REQUIRE( &arr[point[0]][point[1]] == &arr(point[0], point[1]) );
	BOOST_REQUIRE( &arr(point[0], point[1]) == &arr.apply(point) );

#if not defined(__circle_build__)
	BOOST_REQUIRE( &arr[point[0]][point[1]] == &std::apply(arr, point) );
	BOOST_REQUIRE( &arr[point[0]][point[1]] == &     apply(arr, point) );
#endif
}

BOOST_AUTO_TEST_CASE(multi_tests_extension_with_tuple) {
	{
		multi::array<double, 2>::extensions_type ext = {3, 4};
		multi::array<double, 2> arr(ext, 44.);
		BOOST_REQUIRE( size(arr) == 3 );
	}
	{
		auto const [en, em] = std::make_tuple(3, 4);
		multi::array<double, 2> arr({en, em}, 44.);
		BOOST_REQUIRE( size(arr) == 3 );
	}
	{
		auto arr = std::apply([](auto const&... szs) {return multi::array<double, 2>({szs...}, 55.);}, std::make_tuple(3, 4));
		BOOST_REQUIRE( size(arr) == 3 );
		BOOST_REQUIRE( std::get<0>(sizes(arr)) == 3 );
		BOOST_REQUIRE( std::get<1>(sizes(arr)) == 4 );
	}
}

BOOST_AUTO_TEST_CASE(multi_test_constness_reference) {
	multi::array<double, 2> const carr({10, 10}, 99.);

	BOOST_REQUIRE( size( carr(1, {0, 3}) ) == 3 );

	BOOST_REQUIRE( carr(1, {0, 3})[1] == 99. );
	static_assert( decltype( carr({0, 3}, 1) )::rank_v == 1 , "!");
	BOOST_REQUIRE( size(carr.sliced(0, 3)) == 3 );

	BOOST_REQUIRE( carr.range({0, 3}).rotated()[1].unrotated().size() == 3 );

	BOOST_REQUIRE( carr({0, 3}, {0, 3})[1][1] == 99. );

	static_assert(not std::is_assignable_v<decltype(carr(1, {0, 3})[1]), double>, "!");

//  none of these lines should compile because m is read-only
//  m(1, {0, 3})[1] = 88.;
//  m({0, 3}, 1)[1] = 77.;
//  m({0, 3}, {0, 3})[1][1] = 66.;
}

#if 1

BOOST_AUTO_TEST_CASE(multi_test_stencil) {
	multi::array<std::string, 2> arr =
		{{"a", "b", "c", "d", "e"},
		 {"f", "g", "h", "f", "g"},
		 {"h", "i", "j", "k", "l"}}
	;

	BOOST_REQUIRE(      size(arr) == 3                                            );
	BOOST_REQUIRE(           arr.num_elements() == 3*5L                           );
	BOOST_REQUIRE(           arr[1][2] == "h"                                     );

	BOOST_REQUIRE(      size(arr          ({1, 3}, {2, 5})) == 2                  );
	BOOST_REQUIRE( extension(arr          ({1, 3}, {2, 5})).start() == 0          );
	BOOST_REQUIRE(           arr          ({1, 3}, {2, 5}).num_elements() == 2*3L );
	BOOST_REQUIRE(           arr          ({1, 3}, {2, 5}).num_elements() == 2*3L );
	BOOST_REQUIRE(           arr          ({1, 3}, {2, 5})[0][0] == "h"           );
	BOOST_REQUIRE(          &arr          ({1, 3}, {2, 5})[0][0] == &arr[1][2]      );

	BOOST_REQUIRE(      size(arr.stenciled({1, 3}, {2, 5})) == 2                  );
	BOOST_REQUIRE( extension(arr.stenciled({1, 3}, {2, 5})).start() == 1          );
	BOOST_REQUIRE(           arr.stenciled({1, 3}, {2, 5}).num_elements() == 2*3L );
	BOOST_REQUIRE(           arr.stenciled({1, 3}, {2, 5}) [1][2] == "h"          );
	BOOST_REQUIRE(          &arr.stenciled({1, 3}, {2, 5}) [1][2] == &arr[1][2]     );

	BOOST_REQUIRE(  arr().elements().size() == arr.num_elements() );

	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements()[0] == &arr(1, 2) );
	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements()[arr({1, 3}, {2, 5}).elements().size() - 1] == &arr(2, 4) );

	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements().front() == &arr(1, 2) );
	BOOST_REQUIRE( &arr({1, 3}, {2, 5}).elements().back()  == &arr(2, 4) );
}

BOOST_AUTO_TEST_CASE(multi_test_elements_1D) {
	multi::array<double, 1> arr = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
	BOOST_REQUIRE( arr.size() == 10 );

	BOOST_REQUIRE(  arr.elements().size() == 10 );
	BOOST_REQUIRE( &arr.elements()[0] == &arr[0] );
	BOOST_REQUIRE( &arr.elements()[9] == &arr[9] );

	BOOST_REQUIRE(      arr.elements().begin() <  arr.elements().end()     );
	BOOST_REQUIRE(      arr.elements().end()   >  arr.elements().begin()   );
	BOOST_REQUIRE(      arr.elements().begin() != arr.elements().end()     );
	BOOST_REQUIRE( not( arr.elements().begin() == arr.elements().end()   ) );

	BOOST_REQUIRE(  arr().elements().begin() <  arr().elements().end() );
	BOOST_REQUIRE(  arr().elements().begin() == arr().elements().begin() );

	BOOST_REQUIRE( arr().elements().begin() <  arr().elements().end() or arr().elements().begin() == arr().elements().end() );
	BOOST_REQUIRE( arr().elements().begin() <= arr().elements().end() );

	BOOST_REQUIRE(  arr().elements().end()  >  arr().elements().begin() );
	BOOST_REQUIRE(  arr().elements().end()  >= arr().elements().begin() );

	arr.elements() = {9., 8., 7., 6., 5., 4., 3., 2., 1., 0.};
	BOOST_REQUIRE( arr[2] == 7. );
	BOOST_REQUIRE( arr.elements()[2] == 7. );
}

BOOST_AUTO_TEST_CASE(multi_test_elements_1D_as_range) {
	multi::array<double, 1> arr = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
	BOOST_REQUIRE( arr.size() == 10 );

	arr().elements() = {9., 8., 7., 6., 5., 4., 3., 2., 1., 0.};
	BOOST_REQUIRE( arr[2] == 7. );
	BOOST_REQUIRE( arr.elements()[2] == 7. );
}

BOOST_AUTO_TEST_CASE(elements_from_init_list_2D) {
	multi::array<double, 2> arr({3, 2});
	arr().elements() = {1., 2., 3., 4., 5., 6.};
	BOOST_REQUIRE(arr[1][0] == 3.);

	arr.elements() = {10., 20., 30., 40., 50., 60.};
	BOOST_REQUIRE(arr[1][0] == 30.);
}

BOOST_AUTO_TEST_CASE(front_back_2D) {
	multi::array<double, 2> arr({3, 4});
	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 0.);

	BOOST_REQUIRE(  arr.front()[2] ==  arr[0][2] );
	BOOST_REQUIRE( &arr.front()[2] == &arr[0][2] );

	BOOST_REQUIRE(  arr.back ()[2] ==  arr[2][2] );
	BOOST_REQUIRE( &arr.back ()[2] == &arr[2][2] );
}

BOOST_AUTO_TEST_CASE(front_back_1D) {
	multi::array<double, 1> arr({30}, double{});
	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 0.);

	BOOST_REQUIRE(  arr.front() ==  arr[ 0] );
	BOOST_REQUIRE( &arr.front() == &arr[ 0] );

	BOOST_REQUIRE(  arr.back () ==  arr[29] );
	BOOST_REQUIRE( &arr.back () == &arr[29] );
}

BOOST_AUTO_TEST_CASE(elements_rvalues) {
	using movable_type = std::vector<double>;
	movable_type movable_value(5., 99.);

	multi::array<movable_type, 1> arr = {movable_value, movable_value, movable_value};
	BOOST_REQUIRE( arr.size() == 3 );

	movable_type front = std::move(arr)[0];

	BOOST_REQUIRE( front == movable_value );
	BOOST_REQUIRE( arr[0].empty()           );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes
	BOOST_REQUIRE( arr[1] == movable_value  );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

	std::move(arr)[1] = movable_value;
}

template<class Array1D>
void assign_elements_from_to(Array1D&& arr, std::deque<std::vector<double>>& dest) {
	std::copy(std::forward<Array1D>(arr).begin(), std::forward<Array1D>(arr).end(), std::back_inserter(dest));
}

BOOST_AUTO_TEST_CASE(elements_rvalues_nomove) {
	using movable_type = std::vector<double>;
	movable_type movable_value(5., 99.);

	multi::array<movable_type, 1> arr = {movable_value, movable_value, movable_value};
	BOOST_REQUIRE( arr.size() == 3 );

	std::deque<std::vector<double>> q1;

	assign_elements_from_to(arr, q1);

	BOOST_REQUIRE( arr[0] == movable_value );

	std::deque<std::vector<double>> q2;

	assign_elements_from_to(std::move(arr), q2);

	BOOST_REQUIRE( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

	BOOST_REQUIRE( q1 == q2 );
}

BOOST_AUTO_TEST_CASE(elements_rvalues_assignment) {
	std::vector<double> vec = {1., 2., 3.};
	std::move(vec) = std::vector<double>{3., 4., 5.};
	std::move(vec)[1] = 99.;  // it compiles  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes
//  std::move(v[1]) = 99.;  // does not compile

//  double a = 5.;
//	std::move(a) = 9.;  // does not compile
//  BOOST_REQUIRE( a == 9. );

	multi::array<double, 1> arr1 = {1., 2., 3.};
	multi::array<double, 1> arr2 = {1., 2., 3.};
	std::move(arr1) = arr2;  // this compiles TODO(correaa) should it?
}

#endif

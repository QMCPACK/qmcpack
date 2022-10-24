// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi comparisons"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(comparison_complex) {
using complex = std::complex<double>;
{
	multi::array<double , 1> arr = {1., 2., 3.};
	multi::array<complex, 1> arr2 = {1., 2., 3.};
	BOOST_REQUIRE( arr[1] == arr2[1] );
	BOOST_REQUIRE( arr == arr2 ); BOOST_REQUIRE( not (arr != arr2) );
	BOOST_REQUIRE( arr2 == arr ); BOOST_REQUIRE( not (arr2 != arr) );
}
{
	multi::array<double , 2> const arr  = {{1., 2., 3.}, {4., 5., 6.}};
	multi::array<complex, 2> const arr2 = {{1., 2., 3.}, {4., 5., 6.}};
	BOOST_REQUIRE( arr[1][1] == arr2[1][1] );
	BOOST_REQUIRE( arr == arr2 ); BOOST_REQUIRE( not (arr != arr2) );
	BOOST_REQUIRE( arr2 == arr ); BOOST_REQUIRE( not (arr2 != arr) );
	BOOST_REQUIRE( std::equal(arr[1].begin(), arr[1].end(), begin(arr2[1]), end(arr2[1])) );
}
}

BOOST_AUTO_TEST_CASE(multi_comparisons_swap) {
	multi::array<double, 3> arr = {
		{{ 1.2,  1.1}, { 2.4, 1.}},
		{{11.2,  3.0}, {34.4, 4.}},
		{{ 1.2,  1.1}, { 2.4, 1.}}
	};
	BOOST_REQUIRE( arr[0] < arr[1] );

	swap( arr[0], arr[1] );
	BOOST_REQUIRE( arr[1] < arr[0] );

	swap( arr[0], arr[1] );
	BOOST_REQUIRE( arr[0] < arr[1] );
}

BOOST_AUTO_TEST_CASE(comparisons_equality) {
	multi::array<double, 3> arr = {
		{{ 1.2,  1.1}, { 2.4, 1.}},
		{{11.2,  3.0}, {34.4, 4.}},
		{{ 1.2,  1.1}, { 2.4, 1.}}
	};

	multi::array_ref <double, 3>  ref(arr.data_elements(), extensions(arr));
	multi::array_cref<double, 3> cref(data_elements(arr) , extensions(arr));

	BOOST_REQUIRE( arr ==  arr ); BOOST_REQUIRE( not (arr !=  arr) );
	BOOST_REQUIRE( ref ==  arr ); BOOST_REQUIRE( not (ref !=  arr) );
	BOOST_REQUIRE( ref == cref ); BOOST_REQUIRE( not (ref != cref) );

	BOOST_REQUIRE( arr[0] ==  arr[2] );
	BOOST_REQUIRE( ref[0] ==  arr[2] );
	BOOST_REQUIRE( ref[0] == cref[2] );

	BOOST_REQUIRE( not ( arr[0] != arr[2]) );
	BOOST_REQUIRE( not ( ref[0] != ref[2]) );

	BOOST_REQUIRE( not ( arr[0] != arr[2]) );
	BOOST_REQUIRE( not ( ref[0] != ref[2]) );
}

BOOST_AUTO_TEST_CASE(comparisons_ordering) {
	multi::array<double, 3> arr = {
		{{ 1.2,  1.1}, { 2.4, 1.}},
		{{11.2,  3.0}, {34.4, 4.}},
		{{ 1.2,  1.1}, { 2.4, 1.}}
	};

	multi::array_ref <double, 3>  ref(arr.data_elements(), extensions(arr));
	multi::array_cref<double, 3> cref(data_elements(arr) , extensions(arr));

	BOOST_REQUIRE(  arr[0]    <=  arr[1] );
	BOOST_REQUIRE(  ref[0]    <=  arr[1] );
	BOOST_REQUIRE( cref[0]    <= cref[1] );

	BOOST_REQUIRE(  arr[0][0] <= arr[0][1] );
	BOOST_REQUIRE( ref[0][0] <= arr[0][1] );

	BOOST_REQUIRE(  arr[1][0][0] == 11.2 );
	BOOST_REQUIRE(  ref[1][0][0] == 11.2 );
	BOOST_REQUIRE( cref[1][0][0] == 11.2 );

	BOOST_REQUIRE(  arr[0][0][0] == 1.2 );
	BOOST_REQUIRE(  ref[0][0][0] == 1.2 );
	BOOST_REQUIRE( cref[0][0][0] == 1.2 );

	swap(ref[0], ref[1]);

	BOOST_REQUIRE(  begin(arr) <  end(arr) );
	BOOST_REQUIRE( cbegin(arr) < cend(arr) );

	BOOST_REQUIRE(  end(arr) -  begin(arr) == size(arr) );
}

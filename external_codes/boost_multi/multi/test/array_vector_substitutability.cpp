// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi legacy adaptor example"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>
#include<iostream>

namespace multi = boost::multi;

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_1(std::vector<double> const& source, DynamicArray& darr) {
	darr = DynamicArray(source);
}

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_2(std::vector<double> const& source, DynamicArray& darr) {
	darr = DynamicArray(source.begin(), source.end());
}

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_3(std::vector<double> const& source, DynamicArray& darr) {
	darr = std::decay_t<decltype(darr)>(source.begin(), source.end()); // or std::decay_t<decltype(da)>(source.begin(), source.end())
}

template<class It, class DynamicArray>   // e.g. std::vector or multi::array
void resize_copy_4(It first, It last, DynamicArray& darr) {
	darr = DynamicArray(first, last); // or std::decay_t<decltype(da)>(source.begin(), source.end())
}

template<class It, class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_5(It first, It last, DynamicArray& darr) {
	darr.assign(first, last);  // or std::decay_t<decltype(da)>(source.begin(), source.end())
}

// void resize_copy_6   ----> see below test_resize_copy_6

BOOST_AUTO_TEST_CASE(test_resize_copy_1) {
	std::vector<double> const source = {0., 1., 2., 3.};

	std::vector<double>     dest_v = {99., 99.};
	multi::array<double, 1> dest_a = {88., 88.};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_1(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	resize_copy_1(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_2) {
	std::vector<double> const source = {0., 1., 2., 3.};

	std::vector<double>     dest_v = {99., 99.};
	multi::array<double, 1> dest_a = {88., 88.};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_2(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	resize_copy_2(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_3) {
	std::vector<double> const source = {0., 1., 2., 3.};

	std::vector<double>     dest_v = {99., 99.};
	multi::array<double, 1> dest_a = {88., 88.};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_3(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	resize_copy_3(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_4) {
	std::vector<double> const source = {0., 1., 2., 3.};

	std::vector<double>     dest_v = {99., 99.};
	multi::array<double, 1> dest_a = {88., 88.};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_4(source.begin(), source.end(), dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	resize_copy_4(source.begin(), source.end(), dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_5) {
	std::vector<double> const source = {0., 1., 2., 3.};

	std::vector<double>     dest_v = {99., 99.};
	multi::array<double, 1> dest_a = {88., 88.};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_5(source.begin(), source.end(), dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	resize_copy_5(source.begin(), source.end(), dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_6) {
	std::vector<double> const source = {0., 1., 2., 3.};

	std::vector<double>     dest_v = {99., 99.};
	multi::array<double, 1> dest_a = {88., 88.};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	{  // look same code as below
		dest_v = decltype(dest_v)(source);
	}

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	{  // look same code as above
		dest_a = decltype(dest_a)(source);
	}

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

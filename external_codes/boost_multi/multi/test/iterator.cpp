// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi iterators"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<numeric>
#include<vector>

namespace multi = boost::multi;

template<class Array> auto take(Array&& array) -> decltype(array[0]) {return array[0];}

BOOST_AUTO_TEST_CASE(iterator_1d) {
	{
		multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{100}}, 99.);
		BOOST_REQUIRE( size(arr) == 100 );
		BOOST_REQUIRE( begin(arr) < end(arr) );
		BOOST_REQUIRE( end(arr) - begin(arr) == size(arr) );

		multi::array<double, 1>::const_iterator cbarr = cbegin(arr);
		multi::array<double, 1>::iterator barr = begin(arr);
		BOOST_REQUIRE( cbarr == barr );

		multi::array<double, 1>::const_iterator cbarr2 = begin(arr);
		BOOST_REQUIRE( cbarr2 == cbarr );
	}
	{
		multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{100}}, 99.);
		BOOST_REQUIRE( size(arr) == 100 );
		BOOST_REQUIRE( begin(arr) < end(arr) );

		auto arr2 = arr.begin();
		multi::array<double, 1>::const_iterator cbb = arr2;
		BOOST_REQUIRE( cbb == arr2 );
		BOOST_REQUIRE( arr2 == cbb );
	}
}

BOOST_AUTO_TEST_CASE(iterator_2d) {
	{
		multi::array<double, 2> arr({120, 140}, 99.);
		BOOST_REQUIRE(      arr.size() == 120 );
	#if not defined(__circle_build__)  // circle 170 crashes
		BOOST_REQUIRE( size(arr)       == 120 );
	#endif
	#if not defined(__circle_build__)  // circle 170 crashes
		BOOST_REQUIRE( arr.cbegin() < arr.cend() );
	#endif
	#if not defined(__circle_build__)  // circle 170 crashes
		BOOST_REQUIRE( arr.cend() - arr.cbegin() == arr.size() );
	#endif
		using iter = multi::array<double, 2>::iterator;
		static_assert( std::is_same_v< iter::element   , double >, "!");
		static_assert( std::is_same_v< iter::value_type, multi::array<double, 1> >, "!");
		static_assert( std::is_same_v< iter::reference, multi::basic_array<double, 1>>, "!");
		static_assert( std::is_same_v< iter::element_ptr, double*>, "!");

		using citer = multi::array<double, 2>::const_iterator;
		static_assert( std::is_same_v< citer::element   , double >, "!");
		static_assert( std::is_same_v< citer::value_type, multi::array<double, 1> >, "!");
		static_assert( std::is_same_v< citer::reference, multi::basic_array<double, 1, double const*>>, "!");
		static_assert( std::is_same_v< citer::element_ptr, double const* >, "!");
	}
	{
		std::vector<double> vec(10000);
		multi::array_ref<double, 2> arr(vec.data(), {100, 100});
		BOOST_REQUIRE(size(arr) == 100);
		begin(arr)[4][3] = 2.;
	}
}

BOOST_AUTO_TEST_CASE(iterator_interface ) {
	multi::array<double, 3> arr = {
		{
			{ 1.2,  1.1}, { 2.4, 1.}
		},
		{
			{11.2,  3.0}, {34.4, 4.}
		},
		{
			{ 1.2,  1.1}, { 2.4, 1.}
		}
	};

	BOOST_REQUIRE( size(arr)==3 and size(arr[0])==2 and size(arr[0][0])==2);
	BOOST_REQUIRE( arr[0][0][1] == 1.1 );

	BOOST_REQUIRE( begin(arr) < end(arr) );
	BOOST_REQUIRE( cbegin(arr) < cend(arr) );
	BOOST_REQUIRE( begin(arr[0]) < end(arr[0]) );
	BOOST_REQUIRE( begin(arr[0]) < end(arr[0]) );

//  BOOST_REQUIRE(( multi::array<double, 3>::reverse_iterator {A.begin()} == rend(A) ));

//  BOOST_REQUIRE( rbegin(A) < rend(A) );

	BOOST_REQUIRE( end(arr) - begin(arr) == size(arr) );
//  BOOST_REQUIRE( rend(A) - rbegin(A) == size(A) );

	BOOST_REQUIRE( size(*begin(arr)) == 2 );
	BOOST_REQUIRE( size(begin(arr)[1]) == 2 );

	BOOST_REQUIRE( &(arr[1][1].begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
	BOOST_REQUIRE( &arr[0][1][0] == &arr[0][1][0] );
	BOOST_REQUIRE( &((*arr.begin())[1][0]) == &arr[0][1][0] );
	BOOST_REQUIRE( &((*arr.begin()).operator[](1)[0]) == &arr[0][1][0] );
	BOOST_REQUIRE( &(arr.begin()->operator[](1)[0]) == &arr[0][1][0] );
	BOOST_REQUIRE( &(arr.begin()->operator[](1).begin()[0]) == &arr[0][1][0] );  // NOLINT(readability-container-data-pointer) test access
	BOOST_REQUIRE( &((arr.begin()+1)->operator[](1).begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
	BOOST_REQUIRE( &((begin(arr)+1)->operator[](1).begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
	BOOST_REQUIRE( &((cbegin(arr)+1)->operator[](1).begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
}

BOOST_AUTO_TEST_CASE(iterator_semantics) {
	multi::array<double, 3> arr = {
		{{ 1.2,  1.1}, { 2.4, 1.}},
		{{11.2,  3.0}, {34.4, 4.}},
		{{ 1.2,  1.1}, { 2.4, 1.}}
	};

	multi::array<double, 3>::iterator it;
	BOOST_REQUIRE(( multi::array<double, 3>::iterator{} == it ));
	BOOST_REQUIRE(( it == multi::array<double, 3>::iterator{} ));

	it = begin(arr);
	BOOST_REQUIRE( it == begin(arr) );

	multi::array<double, 3>::iterator it2 = begin(arr);
	BOOST_REQUIRE(it == it2);

	it = end(arr);
	BOOST_REQUIRE(it != it2);
	BOOST_REQUIRE(it > it2);

	multi::array<double, 3>::iterator it3{it};
	BOOST_REQUIRE( it3 == it );

	multi::array<double, 3>::const_iterator cit;
	static_assert( std::is_same<multi::array<double, 3>::iterator::element_ptr, double*>{}, "!");

	cit = it3;
	BOOST_REQUIRE( cit == it3 );
	BOOST_REQUIRE( &arr[0][2][1] == &begin(arr)[0][2][1] );

	static_assert( decltype(begin(arr))::rank_v  == 3 , "!" );
	static_assert( decltype(begin(arr))::rank {} == 3 , "!" );

#if not defined(__circle_build__)  // circle 170 crashes
	auto&& ref = multi::ref(begin(arr), end(arr));
//  BOOST_TEST( arr.base() == ref.base() );  // fails in circle (?)
	BOOST_TEST(  arr[0][2][1] ==  ref[0][2][1] );
	BOOST_TEST( &arr[0][2][1] == &ref[0][2][1] );
	BOOST_TEST( arr.layout().stride() == ref.layout().stride());
	BOOST_TEST( arr.layout().offset() == ref.layout().offset());
	BOOST_TEST( arr.layout().nelems() == ref.layout().nelems());

	BOOST_REQUIRE( arr.num_elements() == ref.num_elements() );
	BOOST_REQUIRE( arr.stride() == ref.stride() );
	BOOST_REQUIRE( arr.layout() == ref.layout() );

	BOOST_REQUIRE( &multi::ref(begin(arr), end(arr)) == &arr );
#endif
}

BOOST_AUTO_TEST_CASE(iterator_arrow_operator) {
	multi::array<std::string, 2> arr = {
		{"00", "01"},
		{"10", "11"},
		{"20", "21"}
	};

	BOOST_REQUIRE( arr[1][0] == "10" );

	BOOST_REQUIRE( std::is_sorted(begin(arr), end(arr)) ); // sorted by rows
	BOOST_REQUIRE( std::is_sorted(begin(arr.rotated()), end(arr.rotated())) ); // sorted by cols

	BOOST_REQUIRE( begin( arr           )->size() == arr[0].size() );
	BOOST_REQUIRE( begin( arr.rotated() )->size() == arr.size() );

	BOOST_REQUIRE( &(begin( arr           )->operator[](1)) == &(arr[0][1]) );
	BOOST_REQUIRE( &(begin( arr.rotated() )->operator[](1)) == &(arr[1][0]) );
}

BOOST_AUTO_TEST_CASE(index_range_iteration) {
	multi::index_range irng{0, 5}; // semiopen interval
	std::ostringstream out;
	std::copy(begin(irng), end(irng), std::ostream_iterator<int>{out, ","});
	BOOST_REQUIRE( out.str() == "0,1,2,3,4," );

	BOOST_REQUIRE( std::accumulate(begin(irng), end(irng), 0) == irng.size()*(irng.size()-1)/2 );

	BOOST_REQUIRE( std::accumulate(begin(irng), end(irng), 0, [](auto&& acc, auto const& elem) {return acc + elem*elem*elem;}) > 0 ); // sum of cubes
}

BOOST_AUTO_TEST_CASE(multi_reverse_iterator_1D) {
	multi::array<double, 1> arr(100, 66.);
	BOOST_REQUIRE( &arr[99] == &*std::make_reverse_iterator(arr.end()) );

	auto rbegin = std::make_reverse_iterator(arr.end());
	rbegin += 100;
	multi::array<double, 1>::iterator begin{rbegin.base()};
	BOOST_REQUIRE( begin  == arr.begin() );
}

BOOST_AUTO_TEST_CASE(multi_reverse_iterator_2D) {
	multi::array<double, 2> arr = {
		{  1.,   2.},
		{ 10.,  20.},
		{100., 200.}
	};
	BOOST_REQUIRE( (*arr.begin())[1] == 2. );
	auto rbegin = std::make_reverse_iterator(arr.end());

	BOOST_TEST( (*rbegin)[1] == 200. );
}

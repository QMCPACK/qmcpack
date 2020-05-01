#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x; exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "Unit Tests for Multi sort"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<algorithm> // stable_sort
#include<vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_ref_stable_sort){

	std::vector<double> v = {1.,2.,3.};
	double d2D[4][5] = {
		{150, 16, 17, 18, 19},
		{ 30,  1,  2,  3,  4}, 
		{100, 11, 12, 13, 14}, 
		{ 50,  6,  7,  8,  9} 
	};
	multi::array_ref<double, 2> d2D_ref(&d2D[0][0], {4, 5});
	
	BOOST_REQUIRE( not std::is_sorted(begin(d2D_ref), end(d2D_ref) ) );
	std::stable_sort( begin(d2D_ref), end(d2D_ref) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D_ref), end(d2D_ref) ) );

	BOOST_REQUIRE( not std::is_sorted( begin(d2D_ref<<1), end(d2D_ref<<1) ) );
	std::stable_sort( begin(d2D_ref<<1), end(d2D_ref<<1) );
	BOOST_REQUIRE( std::is_sorted( begin(d2D_ref<<1), end(d2D_ref<<1) ) );

}


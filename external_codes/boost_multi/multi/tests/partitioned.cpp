#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi partitioned operation"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_partitioned_1d){
	multi::array<double, 1>	A1 = {0, 1, 2, 3, 4, 5};
	auto&& A2_ref = A1.partitioned(2);
	static_assert( std::decay<decltype(A2_ref)>::type::dimensionality == decltype(A1)::dimensionality+1 , "!");
	BOOST_REQUIRE( dimensionality(A2_ref)==dimensionality(A1)+1 );
	BOOST_REQUIRE( size(A2_ref)==2 );
	BOOST_REQUIRE( size(A2_ref[0])==3 );
	BOOST_REQUIRE( &A2_ref[1][0] == &A1[3] );
}

BOOST_AUTO_TEST_CASE(array_partitioned_2d){
	multi::array<double, 2>	A2 = 
		{
			{  0,  1,  2,  3,  4,  5}, 
			{  6,  7,  8,  9, 10, 11}, 

			{ 12, 13, 14, 15, 16, 17}, 
			{ 18, 19, 20, 21, 22, 23}, 
		}
	;
	auto&& A3_ref = A2.partitioned(2);
	BOOST_REQUIRE( dimensionality(A3_ref) == dimensionality(A2)+1 );
	BOOST_REQUIRE( num_elements(A3_ref) == num_elements(A2) );
	BOOST_REQUIRE( size(A3_ref)==2 );
	BOOST_REQUIRE( size(A3_ref[0])==2 );
	BOOST_REQUIRE( size(A3_ref[0][0])==6 );
	BOOST_REQUIRE( &A3_ref[1][1][0] == &A2[3][0] );
}

BOOST_AUTO_TEST_CASE(array_partitioned){
	multi::array<std::string, 2> A2 = 
		{
			{  "s0P0",  "s1P0"},
			{  "s0P1",  "s1P1"},
			{  "s0P2",  "s1P2"},
			{  "s0P3",  "s1P3"},
			{  "s0P4",  "s1P4"},
			{  "s0P5",  "s1P5"},
		}; 

	BOOST_REQUIRE(  size(A2) == 6 );
	BOOST_REQUIRE(( sizes(A2) == decltype(sizes(A2)){6, 2} ));

//	BOOST_REQUIRE( size(*(&A2/3)[0]) == 2 );
	BOOST_REQUIRE( size(A2.partitioned(3)) == 3 );
	BOOST_REQUIRE( dimensionality(A2.partitioned(3)) == 3 );

	BOOST_REQUIRE(( sizes(A2.partitioned(3)) == decltype(sizes(A2.partitioned(3))){3, 2, 2} ));
	
}


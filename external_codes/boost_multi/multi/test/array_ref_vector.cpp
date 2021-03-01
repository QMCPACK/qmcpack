#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi array_ref vector"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<iostream>
#include<cassert>
#include<vector>

#include "../array.hpp"

using std::cout; using std::cerr;
namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_ref_vector){ // TODO deside what to do with generalized pointers
//	std::vector<double> buffer(100);
//	multi::array_ref<double, 2, std::vector<double>::iterator> A(buffer.begin(), {10, 10});
//	A[1][1] = 9;
//	BOOST_REQUIRE( A[1][1] == 9 );
//	BOOST_REQUIRE( buffer[11]==9 );

//	A[2]; // requires operator+ 
//	A[1][1]; // requires operator*
}


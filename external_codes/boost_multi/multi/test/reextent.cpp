#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2019

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reextent"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_reextent){
	multi::array<double, 2> A({2, 3});
	BOOST_REQUIRE( num_elements(A) == 6 );

	A[1][2] = 6.;
	BOOST_REQUIRE( A[1][2] == 6. );

	multi::array<double, 2> C({2, 3}); 
	BOOST_REQUIRE(size(C) == 2);
	BOOST_REQUIRE(size(C[0]) == 3);

	A.reextent({5, 4}, 99.); 
	BOOST_REQUIRE( num_elements(A)== 20 );
	BOOST_REQUIRE( A[1][2] == 6. );  // reextent preserves values when it can...
	BOOST_REQUIRE( A[4][3] == 99. ); // ...and gives selected value to the rest

	A = multi::array<double, 2>(extensions(A), 123.); // this is not inefficient, it moves
	BOOST_REQUIRE( A[1][2] == 123. );

	clear(A); // A.clear();
	BOOST_REQUIRE( num_elements(A) == 0 );
	BOOST_REQUIRE( size(A) == 0 );

	A.reextent({5, 4}, 66.);
	BOOST_REQUIRE( A[4][3] == 66. );
}


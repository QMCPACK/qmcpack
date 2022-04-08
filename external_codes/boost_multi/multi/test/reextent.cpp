#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2019

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reextent"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_reextent) {
	multi::array<double, 2> A({2, 3});
	BOOST_REQUIRE( num_elements(A) == 6 );

	A[1][2] = 6.;
	BOOST_REQUIRE( A[1][2] == 6. );

	multi::array<double, 2> C({2, 3});
	BOOST_REQUIRE(size(C) == 2);
	BOOST_REQUIRE(size(C[0]) == 3);

	A.reextent({5, 4}, 99.);
	BOOST_REQUIRE( num_elements(A)== 20 );
	BOOST_TEST_REQUIRE( A[1][2] == 6. );  // reextent preserves values when it can...
	BOOST_REQUIRE( A[4][3] == 99. ); // ...and gives selected value to the rest

	A = multi::array<double, 2>(extensions(A), 123.); // this is not inefficient, it moves
	BOOST_REQUIRE( A[1][2] == 123. );

	clear(A); // A.clear();
	BOOST_REQUIRE( num_elements(A) == 0 );
	BOOST_REQUIRE( size(A) == 0 );

	A.reextent({5, 4}, 66.);
	BOOST_REQUIRE( A[4][3] == 66. );
}

BOOST_AUTO_TEST_CASE(array_reextent_1d) {
	multi::array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}}, 4.);
	BOOST_REQUIRE( size(A) == 10 );
	BOOST_REQUIRE( A[9] == 4. );

	A.reextent(multi::extensions_t<1>{multi::iextension{20}});
	BOOST_REQUIRE( size(A) == 20 );
	BOOST_REQUIRE( A[9] == 4. );
//	BOOST_REQUIRE( A[19] == 0. ); // impossible to know by sometimes 0.

	A.reextent(std::tuple<int>(22) );
	BOOST_REQUIRE( size(A) == 22 );
	BOOST_REQUIRE( A[9] == 4. );


}

BOOST_AUTO_TEST_CASE(array_reextent_0D) {
	multi::array<double, 0> A({}, 4.);
//	A.reextent(A.extensions()); // TODO(correaa) : fix unused for D = 0
	BOOST_REQUIRE( *A.data_elements() == 4. );
}

BOOST_AUTO_TEST_CASE(array_reextent_1d_with_initialization) {
	multi::array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}}, 4.);
	BOOST_REQUIRE( size(A) == 10 );
	BOOST_REQUIRE( A[9] == 4. );

	A.reextent(multi::extensions_t<1>{multi::iextension{20}}, 8.);
	BOOST_REQUIRE( size(A) == 20 );
	BOOST_REQUIRE( A[9] == 4. );
	BOOST_REQUIRE( A[19] == 8. );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d) {
	multi::array<double, 2> A({10, 20}, 4.);
	BOOST_REQUIRE( A[1][2] == 4. );

	A.clear();
	BOOST_REQUIRE( num_elements(A) == 0 );
	BOOST_REQUIRE( size(A) == 0 );

	A.reextent({20, 30}, 9.);
	BOOST_REQUIRE( A[1][2] = 9. );
	BOOST_REQUIRE( A[11][22] = 9. );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d_array) {
	multi::array<double, 2> A({10, 20}, 4.);
	BOOST_REQUIRE( A[1][2] == 4. );

	A.clear();
	BOOST_REQUIRE( num_elements(A) == 0 );
	BOOST_REQUIRE( size(A) == 0 );
}


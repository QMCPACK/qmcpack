// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reversed"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_reversed_3d) {
	multi::array<double, 3> A({30, 40, 50});

	BOOST_TEST_REQUIRE( A.reversed().size() == 50 );

	BOOST_REQUIRE( & A.reversed()[3][5][7] == &A[7][5][3] );
}

template<class Array>
auto flatted_last(Array&& arr) {
	return reversed(flatted(transposed(reversed(std::forward<Array>(arr)))));
}

template<class Array>
auto partitioned_last(Array&& arr, multi::size_type n) {
	return reversed(transposed(partitioned(reversed(std::forward<Array>(arr)), n)));
}

BOOST_AUTO_TEST_CASE(multi_reversed_4d) {
	multi::array<double, 4> A({13, 5, 7, 11});

	BOOST_TEST_REQUIRE( A.reversed().size() == 11 );

	BOOST_REQUIRE( &A.reversed()[1][2][3][4] == &A[4][3][2][1] );

	BOOST_REQUIRE(( sizes(A.reversed().transposed().flatted().reversed()) == decltype(sizes(A.reversed().transposed().flatted().reversed())){13, 5, 77} ));

	BOOST_REQUIRE( &A.reversed().transposed().flatted().reversed()[1][2][ 5] == & A[1][2][0][ 5] );
	BOOST_REQUIRE( &A.reversed().transposed().flatted().reversed()[1][2][10] == & A[1][2][0][10] );
	BOOST_REQUIRE( &A.reversed().transposed().flatted().reversed()[1][2][11] == & A[1][2][1][ 0] );
	BOOST_REQUIRE( &A.reversed().transposed().flatted().reversed()[1][2][12] == & A[1][2][1][ 1] );

	BOOST_REQUIRE( & flatted_last(A)[1][2][12] == & A[1][2][1][1] );
}

BOOST_AUTO_TEST_CASE(multi_reversed_4d_partition_last) {
	multi::array<double, 4> A({11, 5, 7, 12});

	BOOST_TEST_REQUIRE( A.reversed().size() == 12 );

	BOOST_REQUIRE( & A.reversed()[1][2][3][4] == &A[4][3][2][1] );

	BOOST_REQUIRE( & A.reversed().partitioned(3).transposed().reversed()[1][2][3][0][1] == & A[1][2][3][1] );
	BOOST_REQUIRE( & A.reversed().partitioned(3).transposed().reversed()[1][2][3][1][0] == & A[1][2][3][4] );
	BOOST_REQUIRE( & A.reversed().partitioned(3).transposed().reversed()[1][2][3][1][1] == & A[1][2][3][5] );

	BOOST_REQUIRE( & partitioned_last(A, 3)[1][2][3][1][1] == & A[1][2][3][5] );
}


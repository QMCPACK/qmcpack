#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi move"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_swap) {
	multi::array<double, 2> A({3,  5}, 99.);
	multi::array<double, 2> B({7, 11}, 88.);
	swap(A, B);
	BOOST_REQUIRE( size(A) == 7 );
	BOOST_REQUIRE( A[1][2] == 88. );
	BOOST_REQUIRE( B[1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_std_swap) {
	multi::array<double, 2> A({3,  5}, 99.);
	multi::array<double, 2> B({7, 11}, 88.);
	using std::swap;
	swap(A, B);
	BOOST_REQUIRE( size(A) == 7 );
	BOOST_REQUIRE( A[1][2] == 88. );
	BOOST_REQUIRE( B[1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_clear) {
	multi::array<double, 2> A({10, 10}, 99.);
	A.clear();
	BOOST_REQUIRE(A.is_empty());
	A.reextent({20, 20}, 99.);
	BOOST_REQUIRE(not A.is_empty());
	clear(A).reextent({30, 30}, 88.);
	BOOST_REQUIRE(A[15][15] == 88.);
}

BOOST_AUTO_TEST_CASE(multi_array_move) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	multi::array<double, 2> B(std::move(Av[0]), std::allocator<double>{});

	BOOST_REQUIRE( is_empty(Av[0]) );
	BOOST_REQUIRE( size(B) == 4 );
	BOOST_REQUIRE( B[1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv; Bv.reserve(Av.size());

	std::move( begin(Av), end(Av), std::back_inserter(Bv) );

	BOOST_REQUIRE( size(Bv) == size(Av) );
	BOOST_REQUIRE( is_empty(Av[4]) );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_reserve) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv; Bv.reserve(Av.size());

//	for(auto& v: Av) Bv.emplace_back(std::move(v), std::allocator<double>{}); // segfaults nvcc 11.0 but not nvcc 11.1
	std::move(begin(Av), end(Av), std::back_inserter(Bv));

	BOOST_REQUIRE( size(Bv) == size(Av) );
	BOOST_REQUIRE( is_empty(Av[4]) );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_move) {
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv = std::move(Av);

	Av.clear();
	BOOST_REQUIRE( size(Av) == 0 );
	BOOST_REQUIRE( size(Bv) == 10 );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}


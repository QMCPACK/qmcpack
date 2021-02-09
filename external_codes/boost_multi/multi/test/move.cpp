#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi move"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<vector>
#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_array_move){
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	multi::array<double, 2> B(std::move(Av[0]), std::allocator<double>{});

	BOOST_REQUIRE( is_empty(Av[0]) );
	BOOST_REQUIRE( size(B) == 4 );
	BOOST_REQUIRE( B[1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector){
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv;

//	for(auto& v: Av) Bv.emplace_back(std::move(v), std::allocator<double>{}); // segfaults nvcc 11.0 but not nvcc 11.1
	for(auto& v: Av) Bv.emplace_back(std::move(v));

	BOOST_REQUIRE( size(Bv) == size(Av) );
	BOOST_REQUIRE( is_empty(Av[4]) );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_reserve){
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv;
	Bv.reserve(Av.size());

//	for(auto& v: Av) Bv.emplace_back(std::move(v), std::allocator<double>{}); // segfaults nvcc 11.0 but not nvcc 11.1
	for(auto& v: Av) Bv.emplace_back(std::move(v));

	BOOST_REQUIRE( size(Bv) == size(Av) );
	BOOST_REQUIRE( is_empty(Av[4]) );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_move){
	std::vector<multi::array<double, 2> > Av(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > Bv = std::move(Av);

	BOOST_REQUIRE( size(Av) == 0 );
	BOOST_REQUIRE( size(Bv) == 10 );
	BOOST_REQUIRE( size(Bv[5]) == 4 );
	BOOST_REQUIRE( Bv[5][1][2] == 99. );
}


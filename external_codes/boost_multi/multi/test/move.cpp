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

BOOST_AUTO_TEST_CASE(multi_array_move_into_vector){
	std::vector<multi::array<double, 2> > A(10, multi::array<double, 2>({4, 5}, 99.));
	std::vector<multi::array<double, 2> > B;
	B.reserve(10);

	for(auto& v: A){
	//	B.emplace_back(multi::array<double, 2>(std::move(v), std::allocator<double>{}));
		B.emplace_back(std::move(v), std::allocator<double>{});
	//	B.emplace_back(multi::array<double, 2>(std::move(v)));
	//	B.emplace_back(std::move(v));//array<...>{std::move(v), some_communicator} );
	}

	BOOST_REQUIRE( size(B) == size(A) );
	BOOST_REQUIRE( is_empty(A[4]) );
	BOOST_REQUIRE( size(B[5]) == 4 );
	BOOST_REQUIRE( B[5][1][2] == 99. );
}


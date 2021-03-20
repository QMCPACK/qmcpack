#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi allocators"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

//#include "../../multi/memory/stack.hpp" //TODO test custom allocator

#include<complex>
#include<iostream>
#include<scoped_allocator>
#include<vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(std_vector_of_arrays){
	std::vector<multi::array<double, 2>> va;
	std::transform(
		begin(multi::iextension(3)), end(multi::iextension(3)),
		std::back_inserter(va),
		[](auto i){return multi::array<double, 2>({i, i}, static_cast<double>(i));}
	);

	BOOST_REQUIRE( size(va[0]) == 0 );
	BOOST_REQUIRE( size(va[1]) == 1 );
	BOOST_REQUIRE( size(va[2]) == 2 );
	BOOST_REQUIRE( va[1] [0][0] == 1 );
	BOOST_REQUIRE( va[2] [0][0] == 2 );

	std::vector<multi::array<double, 2>> wa = {
		multi::array<double, 2>({0, 0}, 0.),
		multi::array<double, 2>({1, 1}, 1.),
		multi::array<double, 2>({2, 2}, 2.),
	};
	BOOST_REQUIRE( size(va) == size(wa) );
	BOOST_REQUIRE( va == wa );

	std::vector<multi::array<double, 2>> ua(3);
	std::transform(
		begin(multi::iextension(ua.size())), end(multi::iextension(ua.size())), begin(ua), 
		[](auto i){return multi::array<double, 2>({i, i}, static_cast<double>(i));}
	);
	BOOST_REQUIRE( ua == va );
}

BOOST_AUTO_TEST_CASE(array1d_of_arrays2d){
	multi::array<multi::array<double, 2>, 1> A(10, multi::array<double, 2>{});
	BOOST_REQUIRE( size(A) == 10 );

	std::transform(
		begin(extension(A)), end(extension(A)), begin(A), 
		[](auto i){return multi::array<double, 2>({i, i}, static_cast<double>(i));}
	);

	BOOST_REQUIRE( size(A[0]) == 0 );
	BOOST_REQUIRE( size(A[1]) == 1 );
	BOOST_REQUIRE( size(A[8]) == 8 );
	BOOST_REQUIRE( A[8][4][4] == 8 );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d){
	multi::array<multi::array<double, 3>, 2> AA({10, 20}, multi::array<double, 3>{});
	for(int i = 0; i != 10; ++i)
		for(int j = 0; j != 20; ++j)
			AA[i][j] = multi::array<double, 3>({i+j, i+j, i+j}, 99.);
//	std::for_each(begin(extension(AA)), end(extension(AA)), 
//		[&AA](auto& row){std::transform(
//			begin(extension(row)), end(extension(row)), begin(row), 
//			[](auto j){return multi::array<double, 3>({i+j, i+j, i+j}, 99.);}
//		);}
//	);
	BOOST_REQUIRE( size(AA[9][19]) == 9 + 19 );
	BOOST_REQUIRE( AA[9][19][1][1][1] == 99. );
}

BOOST_AUTO_TEST_CASE(array_3d_with_hint_int){
	multi::array<double, 2> const A({3, 4});
	multi::array<int, 3> B({3, 4, 5}, A.cbase());

	B[1][2][3] = 4;
	BOOST_REQUIRE( size(B) == 3 );
	BOOST_REQUIRE( B[1][2][3] == 4 );

	multi::array<int, 3> C({3, 4, 5}, 0);
	BOOST_REQUIRE( size(C) == 3 );

	multi::array<int, 3> D({3, 4, 5}, 99);
	BOOST_REQUIRE( size(D) == 3 );
}

BOOST_AUTO_TEST_CASE(array_3d_with_hint_size_t){
	multi::array<double, 2> const A({3, 4});
	multi::array<size_t, 3> B({3, 4, 5}, A.cbase());

	B[1][2][3] = 4;
	BOOST_REQUIRE( size(B) == 3 );
	BOOST_REQUIRE( B[1][2][3] == 4 );

	multi::array<size_t, 3> C({3, 4, 5}, 0UL);
	BOOST_REQUIRE( size(C) == 3 );

	multi::array<size_t, 3> D({3, 4, 5}, 99);
	BOOST_REQUIRE( size(D) == 3 );
}


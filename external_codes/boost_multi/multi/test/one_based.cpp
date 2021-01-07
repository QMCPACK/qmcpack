#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo Correa 2019-2020

#define BOOST_TEST_DYN_LINK 
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi one based"
#include<boost/test/unit_test.hpp>

#include<iostream>

#include "../array.hpp"
//#include "../adaptors/cuda.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(one_based_1D){

	multi::array<double, 1> Af({{1, 1 + 10}}, 0.);
	Af[1] = 1.;
	Af[2] = 2.;
	Af[3] = 3.;

	BOOST_REQUIRE( Af[1] = 1. );
	BOOST_REQUIRE( *Af.data_elements() == 1. );
	BOOST_REQUIRE( size(Af) == 10 );
	BOOST_REQUIRE( extension(Af).start() == 1 );
	BOOST_REQUIRE( extension(Af).finish() == 11 );

	auto Af1 = multi::array<double, 1>(10, 0.).reindex(1);

	multi::array<double, 1> B({{0, 10}}, 0.);
	B[0] = 1.;
	B[1] = 2.;
	B[2] = 3.;

	BOOST_REQUIRE( size(B) == 10 );
	BOOST_REQUIRE( B != Af );
	BOOST_REQUIRE( std::equal(begin(Af), end(Af), begin(B)) );

	BOOST_REQUIRE( Af.reindexed(0) == B );
}

BOOST_AUTO_TEST_CASE(one_based_2D){

	multi::array<double, 2> Af({{1, 1 + 10}, {1, 1 + 20}}, 0.);
	Af[1][1] = 1.;
	Af[2][2] = 2.;
	Af[3][3] = 3.;
	Af[10][20] = 99.;

	BOOST_REQUIRE( Af[1][1] = 1. );
	BOOST_REQUIRE( Af[10][20] == 99. );
	BOOST_REQUIRE( *Af.data_elements() == 1. );
	BOOST_REQUIRE( Af.data_elements()[Af.num_elements()-1] == 99. );
	BOOST_REQUIRE( size(Af) == 10 );
	BOOST_REQUIRE( extension(Af).start()  ==  1 );
	BOOST_REQUIRE( extension(Af).finish() == 11 );

	auto Af1 = multi::array<double, 2>({10, 10}, 0.).reindex(1, 1);

	multi::array<double, 2> B({{0, 10}, {0, 20}}, 0.);
	B[0][0] = 1.;
	B[1][1] = 2.;
	B[2][2] = 3.;
	B[9][19] = 99.;
	
	BOOST_REQUIRE( size(B) == 10 );
	BOOST_REQUIRE( B != Af );
	BOOST_REQUIRE( std::equal(begin(Af.reindexed(0, 0)), end(Af.reindexed(0, 0)), begin(B)) );
	BOOST_REQUIRE( std::equal(begin(Af), end(Af), begin(B.reindexed(1, 1))) );
	BOOST_REQUIRE( std::equal(begin(Af), end(Af), begin(B.reindexed(0, 1))) );
	
	BOOST_REQUIRE( Af.reindexed(0, 0) == B );

	B = Af;
	BOOST_REQUIRE( B[1][1] = 1. );
	BOOST_REQUIRE( B[10][20] == 99. );
	BOOST_REQUIRE( B == Af );
}

BOOST_AUTO_TEST_CASE(one_base_2D_ref){
	
	double A[3][5] = {
		{ 1.,  2.,  3.,  4.,  5.},
		{ 6.,  7.,  8.,  9., 10.},
		{11., 12., 13., 14., 15.}
	};
	
	multi::array_ref<double, 2> const& Ar = *multi::array_ptr<double, 2>(&A[0][0], {3, 5});
	BOOST_REQUIRE( &Ar[1][3] == &A[1][3] );

	multi::array_ref<double, 2> const& Ar2 = *multi::array_ptr<double, 2>(&A[0][0], {{1, 1+3}, {1, 1+5}});
	BOOST_REQUIRE( sizes(Ar) == sizes(Ar2) );
	BOOST_REQUIRE( &Ar2[1][1] == &A[0][0] );
	BOOST_REQUIRE( &Ar2[2][4] == &A[1][3] );
	
	BOOST_REQUIRE( Ar2 != Ar );
	BOOST_REQUIRE( extensions(Ar2.reindexed(0, 0)) == extensions(Ar) );
	BOOST_REQUIRE( Ar2.reindexed(0, 0) == Ar );
	
	static_assert( not std::is_assignable<decltype(Ar2.reindexed(0, 0)[0][0]), double>{}, "!" );

}


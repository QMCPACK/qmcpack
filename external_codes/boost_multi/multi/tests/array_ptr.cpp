#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi array pointer"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

namespace multi = boost::multi;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(multi_array_ptr, *utf::timeout(2)){
	{
		double a[4][5] = {
			{ 0,  1,  2,  3,  4}, 
			{ 5,  6,  7,  8,  9}, 
			{10, 11, 12, 13, 14}, 
			{15, 16, 17, 18, 19}
		};
		double b[4][5];

		multi::array_ptr<double, 2> aP = &a; // = multi::addressof(a);
		BOOST_REQUIRE( aP->extensions() == multi::extensions(a) );
		BOOST_REQUIRE( extensions(*aP) == multi::extensions(a) );
		using multi::extensions;
		BOOST_REQUIRE( extensions(*aP) == extensions(a) );
		BOOST_REQUIRE( &aP->operator[](1)[1] == &a[1][1] );

		multi::array_ptr<double, 2> aP2 = &a;
		BOOST_REQUIRE( aP == aP2 );

		multi::array_ptr<double, 2> bP = &b;
		BOOST_REQUIRE( bP != aP );

		bP = aP;
		BOOST_REQUIRE( aP == bP );
		BOOST_REQUIRE( *aP == *bP );
		BOOST_REQUIRE( aP->operator==(*bP) );

		auto&& aR = *aP;
		BOOST_REQUIRE( &aR[1][1] == &a[1][1] );
		BOOST_REQUIRE( aR == *aP );
		BOOST_REQUIRE( aR.equal(aP->begin()) );
		BOOST_REQUIRE( size(aR) == aP->size() );
	}
	{
		double a[4][5] = {
			{ 0,  1,  2,  3,  4}, 
			{ 5,  6,  7,  8,  9}, 
			{10, 11, 12, 13, 14}, 
			{15, 16, 17, 18, 19}
		};
		
		std::vector<multi::array_ptr<double, 1>> ps;
		ps.emplace_back(&a[0][0], 5);
		ps.emplace_back(&a[2][0], 5);
		ps.emplace_back(&a[3][0], 5);

		BOOST_TEST( &(*ps[2])[4] == &a[3][4] );
		BOOST_TEST( (*ps[2])[4] == 19 );
		BOOST_TEST( ps[2]->operator[](4) == 19 );

	}
	{
		std::vector<double> v1(100, 3.);
		std::vector<double> const v2(100, 4.);
		multi::array_ptr<double, 2> v1P2D(v1.data(), {10, 10});
		multi::array_cptr<double, 2> v2P2D(v2.data(), {10, 10});

		*v1P2D = *v2P2D;
		v1P2D->operator=(*v2P2D);

		BOOST_REQUIRE( v1[8] == 4. );
	}
}


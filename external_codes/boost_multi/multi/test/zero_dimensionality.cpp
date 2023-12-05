// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi zero dimensionality"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(zero_dimensionality_part1) {
	{
		std::vector<double> v1 = {1., 2., 3.};

		multi::array_ref<double, 1> m1(v1.data(), multi::extensions_t<1>{multi::iextension{3}});
		BOOST_REQUIRE( size(m1) == 3 );
		BOOST_REQUIRE( &m1[1] == &v1[1] );
		BOOST_REQUIRE( num_elements(m1) == 3 );

		multi::array_ref<double, 0> m0(v1.data(), {});
//		BOOST_REQUIRE(( &m0 == multi::array_ptr<double, 0>(v1.data(), {}) ));
		BOOST_REQUIRE( data_elements(m0) == v1.data() );
		BOOST_REQUIRE( num_elements(m0) == 1 );

		m0 = 5.1;
		BOOST_REQUIRE( v1[0] == 5.1 );

		double const& doub = std::move(m0);
		BOOST_REQUIRE( doub == 5.1 );
	}
	{
		multi::static_array<double, 0> a0 = multi::static_array<double, 0>{45.};  // TODO(correaa) this might trigger a compiler crash with g++ 7.5 because of operator&() && overloads
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45. );

		a0 = multi::static_array<double, 0>{60.};
		BOOST_REQUIRE( a0 == 60. );
	}
	{
		std::allocator<double> alloc;
		multi::static_array<double, 0> a0(45., alloc);
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45. );

		a0 = multi::static_array<double, 0>{60.};
		BOOST_REQUIRE( a0 == 60. );
	}
}

BOOST_AUTO_TEST_CASE(zero_dimensionality_part2) {
	{
		multi::array<std::complex<double>, 2> arr({1, 2}, std::allocator<std::complex<double>>{});
		BOOST_REQUIRE( size(arr) == 1 );
	}
	{
		double doub = 2.;
		double dd{multi::array_ref<double, 0>(&doub, {})};

		BOOST_REQUIRE( dd == doub );

		multi::array_ptr<double, 1> ap1(&doub, multi::extensions_t<1>{{0, 1}});
		BOOST_REQUIRE( ap1->base() == &doub );
		BOOST_REQUIRE( (*ap1).base() == &doub );

		multi::array_ptr<double, 0> ap0(&doub, {});

		BOOST_REQUIRE(( ap0 == multi::array_ptr<double, 0>(&doub, {}) ));
		BOOST_REQUIRE(( ap0 != multi::array_ptr<double, 0>(&dd, {}) ));
		BOOST_REQUIRE( ap0->base() == &doub );
		BOOST_REQUIRE( (*ap0).base() == &doub );

		multi::array_ptr<double, 0> ap0dd{&dd};
		BOOST_REQUIRE( ap0dd != ap0 );
		BOOST_REQUIRE( *ap0 == *ap0dd );
		double d3 = M_PI;
		BOOST_REQUIRE(( *multi::array_ptr<double, 0>(&d3, {}) == M_PI ));
	}
}


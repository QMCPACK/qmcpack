#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi References"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<vector>
#include<complex>

namespace multi = boost::multi;
using complex = std::complex<double>;

BOOST_AUTO_TEST_CASE(array_cref_no_mention_type){

	std::vector<double> v(20, 1.);// {1., 2.});

	auto&& v2D = *(v.data()/multi::iextensions<2>{4, 5});
	BOOST_REQUIRE( &v[2] == &v2D[0][2] );
	v2D[1][1] = 30.;

}

BOOST_AUTO_TEST_CASE(array_cref){

	static_assert(std::is_same<std::pointer_traits<complex*>::element_type, complex>{}, "!");
	static_assert(std::is_same<std::pointer_traits<complex*>::rebind<complex const>, complex const*>{}, "!");

	std::vector<complex> d(100);
	std::vector<complex> const dc(100);

	multi::array_ref<complex, 2> A2D(d.data(), multi::iextensions<2>{10, 10});
	multi::array_ref<complex, 2, complex*> B2D(d.data(), {10, 10});
	
	BOOST_REQUIRE( &A2D[3][4] == &B2D[3][4] );
	
//	multi::array_ref<std::complex<double>, 2> C2D(dc.data(), {10, 10}); // error double const* -> double*
	multi::array_ref<complex, 2, complex const*> D2D(dc.data(), {10, 10});
	multi::array_cref<complex, 2> E2D(dc.data(), {10, 10});
	multi::array_cref<complex, 2> F2D(d.data(), {10, 10});
//	F2D[3][4] = 4.; // error, not assignable
	A2D[7][8] = 3.;
	BOOST_REQUIRE( F2D[7][8] == 3. );
	BOOST_REQUIRE( &A2D[7][8] == &F2D[7][8] );

#if defined(__cpp_deduction_guides)
	multi::array_ref G2D(dc.data(), {10, 10}); 
	BOOST_REQUIRE( G2D == D2D );
#endif
	auto&& H2D = multi::make_array_ref<2>(dc.data(), {10, 10}); 
	BOOST_REQUIRE(( H2D == *(dc.data()/multi::iextensions<2>{10, 10}) ));

}


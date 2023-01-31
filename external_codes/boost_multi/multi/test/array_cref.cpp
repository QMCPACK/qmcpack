// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi References"
#include<boost/test/unit_test.hpp>

#include "multi/array_ref.hpp"

#include<complex>
#include<vector>

namespace multi = boost::multi;
using complex = std::complex<double>;

BOOST_AUTO_TEST_CASE(array_cref) {
	static_assert( std::is_same<std::pointer_traits<complex*>::element_type, complex>{}, "!");
	static_assert( std::is_same<std::pointer_traits<complex*>::rebind<complex const>, complex const*>{}, "!");

	std::vector<complex>        vec(100, 0.);
	std::vector<complex> const cvec(100);

	multi::array_ref<complex, 2> A2D(vec.data(), multi::extensions_t<2>{10, 10});
	multi::array_ref<complex, 2, complex*> B2D(vec.data(), {10, 10});

	static_assert( multi::array_ref<complex, 2>::rank::value == 2 , "!" );

	BOOST_REQUIRE( &A2D[3][4] == &B2D[3][4] );

	multi::array_ref <complex, 2, complex const*> D2D(cvec.data(), {10, 10});
	multi::array_cref<complex, 2>                 F2D( vec.data(), {10, 10});

	BOOST_REQUIRE( D2D.layout() == F2D.layout() );

	A2D[7][8] = 3.;
	BOOST_REQUIRE(  F2D[7][8] == 3. );
	BOOST_REQUIRE( &A2D[7][8] == &F2D[7][8] );

//#if defined(__cpp_deduction_guides) and not defined(__NVCC__)
//	multi::array_ref G2D(dc.data(), {10, 10});  // TODO(correaa)
//	BOOST_REQUIRE( G2D == D2D );
//#endif
}


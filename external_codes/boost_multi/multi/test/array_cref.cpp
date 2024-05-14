// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi References"  // NOLINT(cppcoreguidelines-macro-usage) title
#include<boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include<complex>
#include<vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_cref) {
	using complex = std::complex<double>;

	static_assert( std::is_same<std::pointer_traits<complex*>::element_type, complex>{}, "!");
	static_assert( std::is_same<std::pointer_traits<complex*>::rebind<complex const>, complex const*>{}, "!");

	std::vector<complex>        vec(100, 0.);  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	std::vector<complex> const cvec(100);  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	multi::array_ref<complex, 2> A2D(vec.data(), multi::extensions_t<2>{10, 10});
	multi::array_ref<complex, 2, complex*> B2D(vec.data(), {10, 10});

	static_assert( multi::array_ref<complex, 2>::rank::value == 2 , "!" );

	BOOST_REQUIRE( &A2D[3][4] == &B2D[3][4] );

	multi::array_ref <complex, 2, complex const*> const D2D(cvec.data(), {10, 10});
	multi::array_cref<complex, 2>                 F2D( vec.data(), {10, 10});

	BOOST_REQUIRE( D2D.layout() == F2D.layout() );

	A2D[7][8] = 3.;
	BOOST_REQUIRE(  F2D[7][8] == 3. );
	BOOST_REQUIRE( &A2D[7][8] == &F2D[7][8] );

//  #if defined(__cpp_deduction_guides) and not defined(__NVCC__)
//  multi::array_ref G2D(dc.data(), {10, 10});  // TODO(correaa)
//  BOOST_REQUIRE( G2D == D2D );
//  #endif
}

BOOST_AUTO_TEST_CASE(arrays_1D_from_carray) {
	double a_c_array[] = {1.0, 2.0, 3.0};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy types
	multi::array     <double, 1> an_array_value          (a_c_array);  // ok, it is a copy
	multi::array_cref<double, 1> an_array_const_reference(a_c_array);  // ok, it is read only reference
	multi::array_ref <double, 1> an_array_reference      (a_c_array);  // ok, it is a reference

	BOOST_REQUIRE( an_array_value          .size() == 3 and an_array_value          [1] == 2.0 );
	BOOST_REQUIRE( an_array_const_reference.size() == 3 and an_array_const_reference[1] == 2.0 );
	BOOST_REQUIRE( an_array_reference      .size() == 3 and an_array_reference      [1] == 2.0 );
}

BOOST_AUTO_TEST_CASE(arrays_1D_from_const_carray) {
	double const a_c_array[] = {1.0, 2.0, 3.0};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy types
	multi::array     <double, 1> an_array_value          (a_c_array);  // ok, it is a copy
	multi::array_cref<double, 1> an_array_const_reference(a_c_array);  // ok, it is read only reference
//  multi::array_ref <double, 1> an_array_reference      (a_c_array);  // not ok, c array is const

	BOOST_REQUIRE( an_array_value          .size() == 3 and an_array_value          [1] == 2.0 );
	BOOST_REQUIRE( an_array_const_reference.size() == 3 and an_array_const_reference[1] == 2.0 );
//  BOOST_REQUIRE( an_array_reference      .size() == 3 and an_array_reference      [1] == 2. );
}

BOOST_AUTO_TEST_CASE(arrays_1D_from_explict_init_list) {
	std::initializer_list<double> const il = {1.0, 2.0, 3.0};
	multi::array     <double, 1> an_array_value          (il);  // ok, it is a copy
	multi::array_cref<double, 1> an_array_const_reference(il);  // ok, it is read only
//  multi::array_ref <double, 1> an_array_reference      ({1., 2., 3.});  // not allowed, the init list elems are const

	BOOST_REQUIRE( an_array_value           .size() == 3 and an_array_value          [1] == 2.0 );
	BOOST_REQUIRE( an_array_const_reference.size() == 3 and an_array_const_reference[1] == 2.0 );
//  BOOST_REQUIRE( an_array_reference      .size() == 3 and an_array_reference      [1] == 2. );
}

BOOST_AUTO_TEST_CASE(arrays_1D_from_explict_auto_init_list) {
	auto il = {1.0, 2.0, 3.0};
	multi::array     <double, 1> an_array_value          (il);  // ok, it is a copy
	multi::array_cref<double, 1> an_array_const_reference(il);  // ok, it is read only
//  multi::array_ref <double, 1> an_array_reference      ({1., 2., 3.});  // not allowed, the init list elems are const

	BOOST_REQUIRE( an_array_value           .size() == 3 and an_array_value          [1] == 2.0 );
	BOOST_REQUIRE( an_array_const_reference.size() == 3 and an_array_const_reference[1] == 2.0 );
//  BOOST_REQUIRE( an_array_reference      .size() == 3 and an_array_reference      [1] == 2. );
}

BOOST_AUTO_TEST_CASE(arrays_1D_from_init_list) {
	multi::array     <double, 1> an_array_value          ({1.0, 2.0, 3.0});  // ok, it is a copy
//  multi::array_cref<double, 1> an_array_const_reference({1.0, 2.0, 3.0});  // not ok, constructor disable because memcheck detects use after scope
//  multi::array_ref <double, 1> an_array_reference      ({1., 2., 3.});  // not allowed, the init list elems are const

	BOOST_REQUIRE( an_array_value          .size() == 3 and an_array_value          [1] == 2. );
//  BOOST_REQUIRE( an_array_const_reference.size() == 3 and an_array_const_reference[1] == 2.0 );
//  BOOST_REQUIRE( an_array_reference      .size() == 3 and an_array_reference      [1] == 2. );
}

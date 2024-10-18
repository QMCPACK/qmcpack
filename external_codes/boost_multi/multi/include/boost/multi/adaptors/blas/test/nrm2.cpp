// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/test/unit_test.hpp>

#include <boost/multi/adaptors/blas.hpp>
#include <boost/multi/array.hpp>

#include <boost/multi/adaptors/complex.hpp>

#include <complex>

namespace multi = boost::multi;

using complex = multi::complex<double>;
constexpr complex I{0.0, 1.0};  // NOLINT(readability-identifier-length) imaginary unit

BOOST_AUTO_TEST_CASE(multi_blas_nrm2) {
	namespace blas = multi::blas;

	// NOLINTNEXTLINE(readability-identifier-length) blas conventional name
	multi::array<double, 2> const A = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	BOOST_REQUIRE( blas::nrm2(A[1]) == std::sqrt(blas::dot(A[1], A[1])) );

	{
		multi::array<complex, 1> const x = {1.0 + 1.0 * I, 3.0 + 2.0 * I, 3.0 + 4.0 * I};  // NOLINT(readability-identifier-length) blas conventional name
		BOOST_REQUIRE( blas::dot(x, x) == (1.0 + 1.0*I)*(1.0 + 1.0*I) + (3.0 + 2.0*I)*(3.0 + 2.0*I) + (3.0 + 4.0*I)*(3.0 + 4.0*I) );
		using std::sqrt;
		BOOST_REQUIRE( blas::nrm2(x) == sqrt(norm(1.0 + 1.0*I) + norm(3.0 + 2.0*I) + norm(3.0 + 4.0*I)) );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_real) {
	namespace blas                   = multi::blas;
	multi::array<double, 2> const cA = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};

	double n = NAN;  // NOLINT(readability-identifier-length) BLAS naming
	blas::nrm2(rotated(cA)[1], n);

	// BOOST_REQUIRE( blas::nrm2(rotated(cA)[1], n) ==  std::sqrt( 2.0*2.0 + 6.0*6.0 + 10.0*10.0) );  // TODO(correaa) nrm2 is returning a pointer?
	BOOST_REQUIRE( n == std::sqrt( 2.0*2.0 + 6.0*6.0 + 10.0*10.0) );
	// BOOST_REQUIRE( blas::nrm2(rotated(cA)[1]) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

	// double n2 = blas::nrm2(rotated(cA)[1]);
	// BOOST_REQUIRE( n == n2 );

	// multi::array<double, 1> R(4);
	// blas::nrm2( rotated(cA)[1], R[2]);
	// BOOST_REQUIRE( R[2] ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

	// multi::array<double, 0> R0;
	// blas::nrm2( rotated(cA)[1], R0);
	// BOOST_REQUIRE( R0 ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

	// BOOST_REQUIRE( blas::nrm2(rotated(cA)[1]) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
}

BOOST_AUTO_TEST_CASE(multi_adaptor_blas_nrm2_operators) {
	multi::array<double, 1> const X = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming

	double n = NAN;  // NOLINT(readability-identifier-length) BLAS naming

	multi::blas::nrm2(X, n);
	BOOST_REQUIRE( n == multi::blas::nrm2(X) );
}

// BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case){
//   using complex = std::complex<double>;
//   multi::array<complex, 2> const cA = {
//       {1.,  2.,  3.,  4.},
//       {5.,  6.,  7.,  8.},
//       {9., 10., 11., 12.}
//   };

//  using multi::blas::nrm2;
//  double n;
//  BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//  BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//}

// #if 0
// BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case_thrust){
//   using complex = thrust::complex<double>;
//   multi::array<complex, 2> const cA = {
//       {1.,  2.,  3.,  4.},
//       {5.,  6.,  7.,  8.},
//       {9., 10., 11., 12.}
//   };

//  using multi::blas::nrm2;
//  double n;
//  BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//  BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//}

// BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case_types){
//   boost::mpl::for_each<boost::mpl::list<
//       std   ::complex<double>,
//       thrust::complex<double>//,
//   //  boost::multi::complex<double> // TODO make this work
//   >>([](auto cplx){
//       multi::array<decltype(cplx), 2> const cA = {
//           {1.,  2.,  3.,  4.},
//           {5.,  6.,  7.,  8.},
//           {9., 10., 11., 12.}
//       };

//      using multi::blas::nrm2;
//      double n;
//      BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//      BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//  });
//}
// #endif

// BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex){
//   using complex = std::complex<double>; complex const I{0,1};
//   multi::array<complex, 2> const cA = {
//       {1.,  2. + 1.*I,  3.,  4.},
//       {5.,  6. + 4.*I,  7.,  8.},
//       {9., 10. - 3.*I, 11., 12.}
//   };

//  using multi::blas::nrm2;
//  double n;
//  BOOST_REQUIRE( nrm2(rotated(cA)[1], n)   == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );
//  BOOST_REQUIRE( nrm2(rotated(cA)[1])      == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );

//  using namespace multi::blas::operators;
//  BOOST_TEST_REQUIRE( (rotated(cA)[1]^-1) == 1/std::sqrt(norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1])) , boost::test_tools::tolerance(1e-15) );
//  BOOST_TEST_REQUIRE( (rotated(cA)[1]^2) == norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) , boost::test_tools::tolerance(1e-15) );
//}

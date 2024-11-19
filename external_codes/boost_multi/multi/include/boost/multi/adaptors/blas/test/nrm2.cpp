// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

// IWYU pragma: no_include "boost/multi/adaptors/blas/core.hpp"  // for context
// IWYU pragma: no_include "boost/multi/adaptors/blas/traits.hpp"  // for blas, multi
#include <boost/multi/adaptors/blas/core.hpp>  // for context
#include <boost/multi/adaptors/blas/dot.hpp>   // for dot, dot_ref, operator==
#include <boost/multi/adaptors/blas/nrm2.hpp>  // for nrm2, nrm2_ref

#include <boost/multi/array.hpp>               // for array, layout_t, impli...

#include <cmath>  // for sqrt, NAN
#include <complex>
#include <iostream>

namespace multi = boost::multi;

using complex = std::complex<double>;

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imaginary unit

	BOOST_AUTO_TEST_CASE(multi_blas_nrm2) {
		namespace blas = multi::blas;

		// NOLINTNEXTLINE(readability-identifier-length) blas conventional name
		multi::array<double, 2> const A = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};
		BOOST_TEST( blas::nrm2(A[1]) == std::sqrt(blas::dot(A[1], A[1])) );

		{
			multi::array<complex, 1> const x = {1.0 + 1.0 * I, 3.0 + 2.0 * I, 3.0 + 4.0 * I};  // NOLINT(readability-identifier-length) blas conventional name

			BOOST_TEST( std::abs(+blas::dot(x, x) - ((1.0 + 1.0*I)*(1.0 + 1.0*I) + (3.0 + 2.0*I)*(3.0 + 2.0*I) + (3.0 + 4.0*I)*(3.0 + 4.0*I))) < 1.0e-8 );

			std::cout << "nrm2 "<< blas::nrm2(x) << " " << std::sqrt(norm(1.0 + 1.0*I) + norm(3.0 + 2.0*I) + norm(3.0 + 4.0*I)) << '\n';
			BOOST_TEST( std::abs( blas::nrm2(x) - std::sqrt(norm(1.0 + 1.0*I) + norm(3.0 + 2.0*I) + norm(3.0 + 4.0*I)) ) < 1.0e-8 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_real) {
		namespace blas = multi::blas;

		multi::array<double, 2> const cA = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};

		double n = NAN;  // NOLINT(readability-identifier-length) BLAS naming
		blas::nrm2(cA.rotated()[1], n);

		// BOOST_TEST( blas::nrm2(rotated(cA)[1], n) ==  std::sqrt( 2.0*2.0 + 6.0*6.0 + 10.0*10.0) );  // TODO(correaa) nrm2 is returning a pointer?
		BOOST_TEST( n == std::sqrt( 2.0*2.0 + 6.0*6.0 + 10.0*10.0) );

		// BOOST_TEST( blas::nrm2(rotated(cA)[1]) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

		// double n2 = blas::nrm2(rotated(cA)[1]);
		// BOOST_TEST( n == n2 );

		// multi::array<double, 1> R(4);
		// blas::nrm2( rotated(cA)[1], R[2]);
		// BOOST_TEST( R[2] ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

		// multi::array<double, 0> R0;
		// blas::nrm2( rotated(cA)[1], R0);
		// BOOST_TEST( R0 ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

		// BOOST_TEST( blas::nrm2(rotated(cA)[1]) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
	}

	BOOST_AUTO_TEST_CASE(multi_adaptor_blas_nrm2_operators) {
		multi::array<double, 1> const X = {1.1, 2.1, 3.1, 4.1};  // NOLINT(readability-identifier-length) BLAS naming

		{
			double n = NAN;  // NOLINT(readability-identifier-length) BLAS naming

			multi::blas::nrm2(X, n);
			BOOST_TEST( n == multi::blas::nrm2(X) );
		}
		{
			double n = NAN;  // NOLINT(readability-identifier-length) BLAS naming

			n = multi::blas::nrm2(X);
			BOOST_TEST( n == multi::blas::nrm2(X) );
		}
		{
			double const n = multi::blas::nrm2(X);  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( n == multi::blas::nrm2(X) );
		}
		{
			multi::array<double, 0> res{0.0};
			multi::array<double, 1> const xx = {1.0, 2.0, 3.0};

			// multi::blas::dot(xx, xx, res);
			// multi::blas::nrm2(xx, res);
			multi::blas::context ctx;
			multi::blas::dot_n(&ctx, xx.begin(), xx.size(), xx.begin(), res.base());
			// multi::blas::nrm2_n(&ctx, xx.begin(), xx.size(), res.base());

			// BOOST_TEST( *res.base() == 1.0*1.0 + 2.0*2.0 + 3.0*3.0 );

			// multi::blas::nrm2(xx, res);

			// BOOST_TEST( *res.base() == 1.0*1.0 + 2.0*2.0 + 3.0*3.0 );
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case){
	//   using complex = std::complex<double>;
	//   multi::array<complex, 2> const cA = {
	//       {1.0,  2.0,  3.0,  4.0},
	//       {5.0,  6.0,  7.0,  8.0},
	//       {9.0, 10.0, 11.0, 12.0}
	//   };
	// }

	//  using multi::blas::nrm2;
	//  double n;
	//  BOOST_TEST( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
	//  BOOST_TEST( nrm2(rotated(cA)[1])    == n );
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
	//  BOOST_TEST( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
	//  BOOST_TEST( nrm2(rotated(cA)[1])    == n );
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
	//      BOOST_TEST( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
	//      BOOST_TEST( nrm2(rotated(cA)[1])    == n );
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
	//  BOOST_TEST( nrm2(rotated(cA)[1], n)   == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );
	//  BOOST_TEST( nrm2(rotated(cA)[1])      == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );

	//  using namespace multi::blas::operators;
	//  BOOST_TEST_REQUIRE( (rotated(cA)[1]^-1) == 1/std::sqrt(norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1])) , boost::test_tools::tolerance(1e-15) );
	//  BOOST_TEST_REQUIRE( (rotated(cA)[1]^2) == norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) , boost::test_tools::tolerance(1e-15) );
	//}
	return boost::report_errors();
}

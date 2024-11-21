// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/blas/axpy.hpp>     // for axpy, operator-
#include <boost/multi/adaptors/blas/core.hpp>     // for context
#include <boost/multi/adaptors/blas/numeric.hpp>  // for imag, real
// #include <boost/multi/adaptors/complex.hpp>       // for complex, operator*
#include <boost/multi/array.hpp>                  // for array, implicit_cast

#include <boost/core/lightweight_test.hpp>

#include <complex>  // for complex, operator+

namespace multi = boost::multi;
namespace blas  = multi::blas;

// using complex = multi::complex<double>;  // test internal implementation of complex (partially formed complex)
using complex = std::complex<double>;  // test internal implementation of complex (partially formed complex)

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

void zero_stride() {
	multi::array<double, 1> y = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) blas naming

	// regular arrays, stride != 0
	multi::array<double, 1> const x = {1.0, 1.0, 1.0};  // NOLINT(readability-identifier-length) blas naming

	multi::blas::axpy(1.0, x.begin(), y);

	BOOST_TEST(y[0] == 2.0);
	BOOST_TEST(y[1] == 3.0);
	BOOST_TEST(y[2] == 4.0);

	// broadcasted array of 0D into 1D, stride == 0
	multi::array<double, 0> const x0(1.0);

	// x0.broadcasted() is {..., 1.0, 1.0, ...} and has zero stride

	multi::blas::axpy(1.0, x0.broadcasted().begin(), y);

	BOOST_TEST(y[0] == 3.0);
	BOOST_TEST(y[1] == 4.0);
	BOOST_TEST(y[2] == 5.0);
}

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(multi_blas_axpy_real) {
		multi::array<double, 2> arr = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};

		auto const AC = arr;

		multi::array<double, 1> const b = arr[2];  // NOLINT(readability-identifier-length) BLAS naming

		blas::axpy(2.0, b, arr[1]);  // daxpy
		BOOST_TEST( arr[1][2] == 2.0*b[2] + AC[1][2] );
	}

	BOOST_AUTO_TEST_CASE(blas_axpy_repeat) {
		multi::array<double, 1> a1D = multi::iextension(3);
		BOOST_TEST( a1D[0] == 0.0 );
		BOOST_TEST( a1D[1] == 1.0 );
		BOOST_TEST( a1D[2] == 2.0 );

		multi::array<double, 1> const b1D = {3.0, 3.0, 3.0};

		blas::axpy(1.0, b1D, a1D);
		BOOST_TEST( a1D[0] == 3.0 );
		BOOST_TEST( a1D[1] == 4.0 );
		BOOST_TEST( a1D[2] == 5.0 );

		// BOOST_TEST(( multi::array<double, 0>(3.0).broadcasted().size() != 0 ));

		blas::axpy_n(1.0, multi::array<double, 0>(3.0).broadcasted().begin(), 3, a1D.begin());
		BOOST_TEST( a1D[0] == 6.0 );
		BOOST_TEST( a1D[1] == 7.0 );
		BOOST_TEST( a1D[2] == 8.0 );

		// blas::axpy(1.0, multi::array<double, 0>(3.0).broadcasted(), a1D);
		// BOOST_TEST( a1D[0] == 6.0 );
		// BOOST_TEST( a1D[1] == 7.0 );
		// BOOST_TEST( a1D[2] == 8.0 );

		// blas::axpy(2.0, b, arr[1]);  // daxpy
		// BOOST_TEST( arr[1][2] == 2.0*b[2] + AC[1][2] );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_axpy_double) {
		multi::array<double, 2> const const_arr = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};
		multi::array<double, 2>       arr = const_arr;
		multi::array<double, 1> const b{const_arr[2]};  // NOLINT(readability-identifier-length) conventional name in BLAS

		blas::axpy(2.0, b, arr[1]);  // A[1] = 2*b + A[1], A[1]+= a*A[1]
		BOOST_TEST( arr[1][2] == 2.0*b[2] + const_arr[1][2] );

		auto const I = complex{0, 1};  // NOLINT(readability-identifier-length) imaginary unit

		multi::array<complex, 1> AC = {1.0 + 2.0 * I, 3.0 + 4.0 * I, 4.0 - 8.0 * I};
		multi::array<complex, 1> BC(extensions(AC), complex{0.0, 0.0});

		blas::axpy(+1.0, blas::real(AC), blas::real(BC));
		blas::axpy(-1.0, blas::imag(AC), blas::imag(BC));

		//  BOOST_TEST( BC[2] == std::conj(AC[2]) );
		BOOST_TEST( BC[2] == conj(AC[2]) );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex) {
		multi::array<complex, 2> arr = {
			{{1.0, 0.0},  {2.0, 0.0},  {3.0, 0.0},  {4.0, 0.0}},
			{{5.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0}},
			{{9.0, 0.0}, {10.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}},
		};
		auto const const_arr = arr;

		multi::array<complex, 1> const x = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
		blas::axpy(complex{2.0, 0.0}, x, arr[1]);   // zaxpy (2. is promoted to 2+I*0 internally and automatically)
		BOOST_TEST( arr[1][2] == 2.0*x[2] + const_arr[1][2] );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_as_operator_plus_equal) {
		using complex = std::complex<double>;

		multi::array<complex, 2> arr = {
			{{1.0, 0.0},  {2.0, 0.0},  {3.0, 0.0},  {4.0, 0.0}},
			{{5.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0}},
			{{9.0, 0.0}, {10.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}},
		};
		auto const                     carr = arr;
		multi::array<complex, 1> const y    = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
		arr[1] += blas::axpy(2.0, y);                  // zaxpy (2. is promoted to 2+I*0 internally and automatically)
		BOOST_TEST( arr[1][2] == 2.0*y[2] + carr[1][2] );
	}

	BOOST_AUTO_TEST_CASE(axpy_assignment) {
		multi::array<double, 1> const xx = {1.0, 1.0, 1.0};
		multi::array<double, 1>       yy = {2.0, 2.0, 2.0};

		yy += blas::axpy(3.0, xx);

		BOOST_TEST( yy[0] == 5.0 );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_as_operator_minus_equal) {
		multi::array<complex, 2> arr = {
			{{1.0, 0.0},  {2.0, 0.0},  {3.0, 0.0},  {4.0, 0.0}},
			{{5.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0}},
			{{9.0, 0.0}, {10.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}},
		};
		auto const                     AC = arr;
		multi::array<complex, 1> const x  = arr[2];  // NOLINT(readability-identifier-length) BLAS naming
		arr[1] -= blas::axpy(complex{2.0, 0.0}, x);  // zaxpy (2. is promoted to 2+I*0 internally and automatically)
		BOOST_TEST( arr[1][2] == -2.0*x[2] + AC[1][2] );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_axpy_complex_context) {
		multi::array<complex, 2> arr = {
			{{1.0, 0.0},  {2.0, 0.0},  {3.0, 0.0},  {4.0, 0.0}},
			{{5.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0}},
			{{9.0, 0.0}, {10.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}},
		};
		auto const                     arr_copy = arr;
		multi::array<complex, 1> const arr2     = arr[2];
		blas::context                  ctxt{};
		blas::axpy(&ctxt, complex{2.0, 0.0}, arr2, arr[1]);  // zaxpy (2. is promoted to 2+I*0 internally and automatically)
		BOOST_TEST( arr[1][2] == 2.0*arr2[2] + arr_copy[1][2] );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_axpy_operator_minus) {
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 1> x = {
			{10.0, 0.0},
			{11.0, 0.0},
			{12.0, 0.0},
			{13.0, 0.0},
		};
		multi::array<complex, 1> const y = x;  // NOLINT(readability-identifier-length) BLAS naming

		using blas::operators::operator-;

		BOOST_TEST( (x - y)[0] == complex(0.0, 0.0) );
		BOOST_TEST( (y - x)[0] == complex(0.0, 0.0) );

		using blas::operators::operator+;

		BOOST_TEST( (x - (y+y))[0] == -x[0] );
		BOOST_TEST( ((x+x) - y)[0] == +x[0] );

		multi::array<complex, 2> arr = {
			{{1.0, 0.0}, {2.0, 0.0}},
			{{3.0, 0.0}, {4.0, 0.0}},
		};
		multi::array<complex, 1> const arr2 = {
			{1.0, 0.0},
			{2.0, 0.0},
		};
		BOOST_TEST( (arr[0] - arr2)[0] == complex(0.0, 0.0) );
		BOOST_TEST( (arr[0] - arr2)[1] == complex(0.0, 0.0) );

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 1> X = {
			{10.0, 0.0},
			{11.0, 0.0},
			{12.0, 0.0},
			{13.0, 0.0},
		};

		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::array<complex, 1> const Y = {
			{10.0, 0.0},
			{11.0, 0.0},
			{12.0, 0.0},
			{13.0, 0.0},
		};

		using blas::operators::operator-=;
		X -= Y;
		BOOST_TEST( X[0] == complex(0.0, 0.0) );
	}

	zero_stride();

	return boost::report_errors();
}

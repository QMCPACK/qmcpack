// Copyright 2019-2024 Alfredo A. Correa

#include <boost/test/unit_test.hpp>

#include "../../../array.hpp"

#include "../../blas/copy.hpp"

#include <complex>

namespace multi = boost::multi;
namespace blas  = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_copy_n) {
	multi::array<double, 1> const x = {1.0, 2.0, 3.0, 4.0};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1>       y = {5.0, 6.0, 7.0, 8.0};  // NOLINT(readability-identifier-length) BLAS naming
	blas::copy_n(x.begin(), x.size(), y.begin());
	BOOST_REQUIRE( y == x );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy) {
	multi::array<double, 1> const x = {1.0, 2.0, 3.0, 4.0};  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<double, 1> y = {5.0, 6.0, 7.0, 8.0};  // NOLINT(readability-identifier-length) BLAS naming
		blas::copy(x, y);  // segmentation fault in clang-11
		BOOST_REQUIRE( y == x );
	}
	{
		multi::array<double, 1> y = {5.0, 6.0, 7.0, 8.0};  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE( size(y) == size(x) );
		y() = blas::copy(x);
		BOOST_REQUIRE( y == x );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_real) {
	namespace blas              = multi::blas;
	multi::array<double, 2> arr = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	BOOST_REQUIRE( arr[0][2] ==  3.0 );
	BOOST_REQUIRE( arr[2][2] == 11.0 );

	blas::copy(arr[0], arr[2]);
	BOOST_REQUIRE( arr[0][2] ==  3.0 );
	BOOST_REQUIRE( arr[2][2] ==  3.0 );

	blas::copy(arr[1]({0, size(arr[1])}), arr[2]({0, size(arr[1])}));
	BOOST_REQUIRE( arr[1][3] == 8.0 );
	BOOST_REQUIRE( arr[2][3] == 8.0 );

	multi::array<double, 1> AR3 = blas::copy(rotated(arr)[3]);  // dcopy
	BOOST_REQUIRE( AR3[1] == arr[1][3] );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy_row) {
	multi::array<double, 2> const arr = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0},
	};
	multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{3}});  // NOLINT(readability-identifier-length) BLAS naming
	blas::copy(rotated(arr)[0], y);
	BOOST_REQUIRE( y == rotated(arr)[0] );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_complex) {
	using complex                = std::complex<double>;
	auto const               I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> arr = {
		{1.0 + 3.0 * I,  2.0 + 4.0 * I,  3.0 + 5.0 * I,  4.0 + 6.0 * I},
		{5.0 + 0.0 * I,  6.0 + 0.0 * I,  7.0 + 0.0 * I,  8.0 + 0.0 * I},
		{9.0 + 0.0 * I, 10.0 + 0.0 * I, 11.0 + 0.0 * I, 12.0 + 0.0 * I},
	};
	blas::copy(arr[0], arr[2]);
	BOOST_REQUIRE( arr[0][2] == 3.0 + 5.0*I );
}

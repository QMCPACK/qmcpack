// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS copy"
#include<boost/test/unit_test.hpp>

#include "../../../array.hpp"

#include "../../blas/copy.hpp"

#include<complex>

#include "config.hpp"

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_copy_n) {
	multi::array<double, 1> const x = {1., 2., 3., 4.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> y = {5., 6., 7., 8.};  // NOLINT(readability-identifier-length) BLAS naming
	blas::copy_n(x.begin(), x.size(), y.begin());
	BOOST_REQUIRE( y == x );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy_it) {
	multi::array<double, 1> const x = {1., 2., 3., 4.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> y = {5., 6., 7., 8.};  // NOLINT(readability-identifier-length) BLAS naming
	blas::copy(x.begin(), x.end(), y.begin());
	BOOST_REQUIRE( y == x );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy) {
	multi::array<double, 1> const x = {1., 2., 3., 4.};  // NOLINT(readability-identifier-length) BLAS naming
	{
		multi::array<double, 1> y = {5., 6., 7., 8.};  // NOLINT(readability-identifier-length) BLAS naming
		blas::copy(x, y);  // segmentation fault in clang-11
		BOOST_REQUIRE( y == x );
	}
	{
		multi::array<double, 1> y = {5., 6., 7., 8.};  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE( size(y) == size(x) );
		y = blas::copy(x);
		BOOST_REQUIRE( y == x );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_real) {
	namespace blas = multi::blas;
	multi::array<double, 2> arr = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	BOOST_REQUIRE( arr[0][2] ==  3. );
	BOOST_REQUIRE( arr[2][2] == 11. );

	blas::copy(arr[0], arr[2]);
	BOOST_REQUIRE( arr[0][2] ==  3. );
	BOOST_REQUIRE( arr[2][2] ==  3. );

//	multi::blas::copy(begin(A[1]), end(A[1]), begin(A[2])); // dcopy
	blas::copy( arr[1]({0, size(arr[1])}), arr[2]({0, size(arr[1])}) );
	BOOST_REQUIRE( arr[1][3] == 8. );
	BOOST_REQUIRE( arr[2][3] == 8. );

	multi::array<double, 1> AR3 = blas::copy(rotated(arr)[3]); // dcopy
	BOOST_REQUIRE( AR3[1] == arr[1][3] );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy_row) {
	multi::array<double, 2> const arr = {
		{1., 2., 3.},
		{4., 5., 6.},
		{7., 8., 9.}
	};
	multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{3}});  // NOLINT(readability-identifier-length) BLAS naming
	blas::copy(rotated(arr)[0], y);
	BOOST_REQUIRE( y == rotated(arr)[0] );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_complex) {
	using complex = std::complex<double>; constexpr complex I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> arr = {
		{1. + 3.*I,  2. + 4.*I,  3. + 5.*I,  4. + 6.*I},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	blas::copy(arr[0], arr[2]);
	BOOST_REQUIRE( arr[0][2] == 3. + 5.*I );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy_context) {
	multi::array<double, 1> const x = {1., 2., 3., 4.};  // NOLINT(readability-identifier-length) BLAS naming
	blas::context ctx;
	{
		multi::array<double, 1> y = {5., 6., 7., 8.};  // NOLINT(readability-identifier-length) BLAS naming
		blas::copy(ctx, x, y);
		BOOST_REQUIRE( x == y );
	}
	{
		multi::array<double, 1> y = {5., 6., 7., 8.};  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE( size(y) == size(x) );
		y = blas::copy(ctx, x);
		BOOST_REQUIRE( x == y );
	}
}

#if CUDA_FOUND
#include<thrust/complex.h>

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_copy_thrust) {
	multi::array<thrust::complex<double>, 1> const a(multi::extensions_t<1>{multi::iextension{10}}, thrust::complex<double>{});
	multi::array<thrust::complex<double>, 1> b(multi::extensions_t<1>{multi::iextension{10}});
	blas::copy(a, b);

	BOOST_REQUIRE( a == b );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_text_copy_interop) {
	static_assert( std::is_convertible<std::complex<double>, thrust::complex<double>>{} );
	static_assert( std::is_convertible<thrust::complex<double>, std::complex<double>>{} );
	multi::array<std::complex<double>, 1> a(multi::extensions_t<1>{multi::iextension{10}}, std::complex<double>{});
	multi::array<thrust::complex<double>, 1> b(multi::extensions_t<1>{multi::iextension{10}});
	blas::copy(a, b);

	BOOST_REQUIRE( a == b );
}
#endif

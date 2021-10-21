// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS copy"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../array.hpp"

#include "../../blas/copy.hpp"

#include<complex>

#include "config.hpp"

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_copy_n) {
	multi::array<double, 1> const A = {1., 2., 3., 4.};
	multi::array<double, 1> B = {5., 6., 7., 8.};
	blas::copy_n(A.begin(), A.size(), B.begin());
	BOOST_REQUIRE( B == A );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy_it) {
	multi::array<double, 1> const A = {1., 2., 3., 4.};
	multi::array<double, 1> B = {5., 6., 7., 8.};
	blas::copy(A.begin(), A.end(), B.begin());
	BOOST_REQUIRE( B == A );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy) {
	multi::array<double, 1> const A = {1., 2., 3., 4.};
	{
		multi::array<double, 1> B = {5., 6., 7., 8.};
		blas::copy(A, B); // segmentation fault in clang-11
		BOOST_REQUIRE( B == A );
	}
	 {
		multi::array<double, 1> B = {5., 6., 7., 8.};
		BOOST_REQUIRE( size(B) == size(A) );
		B = blas::copy(A);
		BOOST_REQUIRE( B == A );
	}
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_real) {
	namespace blas = multi::blas;
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	BOOST_REQUIRE( A[0][2] ==  3. );
	BOOST_REQUIRE( A[2][2] == 11. );

	blas::copy(A[0], A[2]);
	BOOST_REQUIRE( A[0][2] ==  3. );
	BOOST_REQUIRE( A[2][2] ==  3. );

//	multi::blas::copy(begin(A[1]), end(A[1]), begin(A[2])); // dcopy
	blas::copy( A[1]({0, size(A[1])}), A[2]({0, size(A[1])}) );
	BOOST_REQUIRE( A[1][3] == 8. );
	BOOST_REQUIRE( A[2][3] == 8. );

	multi::array<double, 1> AR3 = blas::copy(rotated(A)[3]); // dcopy
	BOOST_REQUIRE( AR3[1] == A[1][3] );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy_row) {
	multi::array<double, 2> const A = {
		{1., 2., 3.},
		{4., 5., 6.},
		{7., 8., 9.}
	};
	multi::array<double, 1> B(multi::extensions_t<1>{multi::iextension{3}});
	blas::copy(rotated(A)[0], B);
	BOOST_REQUIRE( B == rotated(A)[0] );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_complex) {
	using complex = std::complex<double>; constexpr complex I{0, 1};
	multi::array<complex, 2> A = {
		{1. + 3.*I,  2. + 4.*I,  3. + 5.*I,  4. + 6.*I},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	blas::copy(A[0], A[2]);
	BOOST_REQUIRE( A[0][2] == 3. + 5.*I );
}

BOOST_AUTO_TEST_CASE(multi_blas_copy_context) {
	multi::array<double, 1> const A = {1., 2., 3., 4.};
	blas::context ctx;
	{
		multi::array<double, 1> B = {5., 6., 7., 8.};
		blas::copy(ctx, A, B);
		BOOST_REQUIRE( A == B );
	}
	 {
		multi::array<double, 1> B = {5., 6., 7., 8.};
		BOOST_REQUIRE( size(B) == size(A) );
		B = blas::copy(ctx, A);
		BOOST_REQUIRE( A == B );
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


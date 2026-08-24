// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2024 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS asum"
// #include <boost/test/unit_test.hpp>

#include "../../../adaptors/cuda.hpp"
#include "../../../array.hpp"
#include "../../blas/asum.hpp"
#include "../../blas/cuda.hpp"

#include "multi/adaptors/complex.hpp"

#include <complex>
#include <numeric>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(const multi_blas_asum_double) {
	multi::array<double, 2> const A = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	using multi::blas::asum;
	BOOST_REQUIRE(asum(A[1]) == std::accumulate(begin(A[1]), end(A[1]), 0.0, [](auto&& a, auto&& b) { return a + std::abs(b); }));
}

BOOST_AUTO_TEST_CASE(const multi_blas_asum_complex) {
	using Z = multi::complex<double>;  // std::complex<double>;

	auto const I = Z{0.0, 1.0};

	multi::array<Z, 2> const A = {
		{1.0 + 2.0 * I,           2.0,            3.0,  4.0},
		{          5.0, 6.0 + 3.0 * I,            7.0,  8.0},
		{          9.0,          10.0, 11.0 + 4.0 * I, 12.0},
	};

	using multi::blas::asum;
	BOOST_REQUIRE(
		asum(A[1]) == std::accumulate(
			begin(A[1]), end(A[1]), 0.0,
			[](auto&& a, auto&& b) { return a + std::abs(real(b)) + std::abs(imag(b)); }
		)
	);
}

BOOST_AUTO_TEST_CASE(const multi_blas_asum_double_cuda) {
	multi::cuda::array<double, 2> const A = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	using multi::blas::asum;
	BOOST_REQUIRE(asum(A[1]) == 26.0 );
}

using complex = multi::complex<double>;
constexpr auto I = complex{0.0, 1.0};

BOOST_AUTO_TEST_CASE(const multi_blas_asum_complex_cuda) {
	namespace blas = multi::blas;

	multi::cuda::array<complex, 2> const A = {
		{1.0 + 2.0 * I,  2.0          ,  3.0          ,  4.0},
		{5.0          ,  6.0 + 3.0 * I,  7.0          ,  8.0},
		{9.0          , 10.0          , 11.0 + 4.0 * I, 12.0},
	};

	BOOST_REQUIRE( blas::asum(A[1]        ) == 29.0 );
	BOOST_REQUIRE( blas::asum(A[1]({0, 4})) == 29.0 );
}

BOOST_AUTO_TEST_CASE(const multi_blas_asum_complex_cuda_mutable) {
	using Z = multi::complex<double>;

	auto const I = Z{0.0, 1.0};

	multi::cuda::array<Z, 2> const A = {
		{1.0 + 2.0 * I,  2.0          ,  3.0          ,  4.0},
		{5.0          ,  6.0 + 3.0 * I,  7.0          ,  8.0},
		{9.0          , 10.0          , 11.0 + 4.0 * I, 12.0},
	};

	using multi::blas::asum;
	BOOST_REQUIRE( asum(A[1]        ) == Z{29.0} );
	BOOST_REQUIRE( asum(A[1]({0, 4})) == Z{29.0} );
}

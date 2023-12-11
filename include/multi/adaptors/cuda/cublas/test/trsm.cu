// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS trsm"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/cuda/cublas.hpp>
#include <multi/adaptors/blas/trsm.hpp>

#include <multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(unit_trsm_multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const A = {
		{ 1.0 + 0.0*I, 3.0 + 1.2*I,  5.0 - 12.0*I},
		{ 0.0 + 0.0*I, 1.0 + 0.0*I,  2.1 +  1.1*I},
		{ 0.0 + 0.0*I, 0.0 + 0.0*I,  1.0 +  0.0*I},
	};
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> B = {
		{1.0 + 1.0*I, 2.0 + 1.0*I, 3.0 + 1.0*I},
		{5.0 + 3.0*I, 9.0 + 3.0*I, 1.0 - 1.0*I},
	};

	multi::thrust::cuda::array<complex, 2> const A_gpu = A;
	multi::thrust::cuda::array<complex, 2> B_gpu = B;

	using multi::blas::trsm;
	using multi::blas::filling;
	using multi::blas::hermitized;

	blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::unit, complex{1.0, 0.0}, A    , blas::H(B));  // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_REQUIRE_CLOSE( B[1][0].real() , -43.439999999999998 , 0.001);
	BOOST_REQUIRE_CLOSE( B[1][0].imag() , -13.000000000000002 , 0.001);

	blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::unit, complex{1.0, 0.0}, A_gpu, blas::H(B_gpu));
	multi::array<complex, 2> B_cpy = B_gpu;

	BOOST_REQUIRE_CLOSE( B_cpy[1][0].real() , -43.439999999999998 , 0.001);
	BOOST_REQUIRE_CLOSE( B_cpy[1][0].imag() , -13.000000000000002 , 0.001);
}

BOOST_AUTO_TEST_CASE(trsm_multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const A = {
		{ 1.0 + 4.0*I, 3.0 + 0.0*I,  4.0 - 10.0*I},
		{ 0.0 + 0.0*I, 7.0 - 3.0*I,  1.0 +  0.0*I},
		{ 0.0 + 0.0*I, 0.0 + 0.0*I,  8.0 -  2.0*I},
	};
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> B = {
		{1.0 + 1.0*I, 2.0 + 1.0*I, 3.0 + 1.0*I},
		{5.0 + 3.0*I, 9.0 + 3.0*I, 1.0 - 1.0*I},
	};

	multi::thrust::cuda::array<complex, 2> const A_gpu = A;
	multi::thrust::cuda::array<complex, 2> B_gpu = B;

	using multi::blas::trsm;
	using multi::blas::filling;
	using multi::blas::hermitized;

	// B = ConjugateTranspose[Inverse[A] . ConjugateTranspose[B]]
	// ConjugateTranspose[B] = Inverse[A] . ConjugateTranspose[B]
	blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::non_unit, complex{1.0, 0.0}, A    , blas::H(B));  // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_REQUIRE_CLOSE( B[1][0].real(), -0.72562939983295538 , 0.001);
	BOOST_REQUIRE_CLOSE( B[1][0].imag(),  0.046772461520104877, 0.001);

	BOOST_REQUIRE_CLOSE( real(blas::H(B)[0][1]), -0.72562939983295538 , 0.001);
	BOOST_REQUIRE_CLOSE( imag(blas::H(B)[0][1]), -0.046772461520104877, 0.001);

	blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::non_unit, complex{1.0, 0.0}, A_gpu, blas::H(B_gpu));
	cudaDeviceSynchronize();

	multi::array<complex, 2> B_cpy = B_gpu;
	BOOST_REQUIRE_CLOSE( B_cpy[1][0].real() , -0.72562939983295538 , 0.001);
	BOOST_REQUIRE_CLOSE( B_cpy[1][0].imag() ,  0.046772461520104877, 0.001);
}

BOOST_AUTO_TEST_CASE(default_param_unit_trsm_multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const A = {
		{ 1.0 + 4.0*I, 3.0 + 0.0*I,  4.0 - 10.0*I},
		{ 0.0 + 0.0*I, 7.0 - 3.0*I,  1.0 +  0.0*I},
		{ 0.0 + 0.0*I, 0.0 + 0.0*I,  8.0 -  2.0*I},
	};
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> B = {
		{1.0 + 1.0*I, 2.0 + 1.0*I, 3.0 + 1.0*I},
		{5.0 + 3.0*I, 9.0 + 3.0*I, 1.0 - 1.0*I},
	};

	multi::thrust::cuda::array<complex, 2> const A_gpu = A;
	multi::thrust::cuda::array<complex, 2> B_gpu = B;

	using multi::blas::trsm;
	using multi::blas::filling;
	using multi::blas::hermitized;

	// B = ConjugateTranspose[Inverse[A] . ConjugateTranspose[B]]
	// ConjugateTranspose[B] = Inverse[A] . ConjugateTranspose[B]
	blas::trsm(blas::side::left, blas::filling::upper, complex{1.0, 0.0}, A    , blas::H(B));  // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_REQUIRE_CLOSE( B[1][0].real(), -0.72562939983295538 , 0.001);
	BOOST_REQUIRE_CLOSE( B[1][0].imag(),  0.046772461520104877, 0.001);

	BOOST_REQUIRE_CLOSE( real(blas::H(B)[0][1]), -0.72562939983295538 , 0.001);
	BOOST_REQUIRE_CLOSE( imag(blas::H(B)[0][1]), -0.046772461520104877, 0.001);

	blas::trsm(blas::side::left, blas::filling::upper, blas::diagonal::non_unit, complex{1.0, 0.0}, A_gpu, blas::H(B_gpu));
	cudaDeviceSynchronize();

	multi::array<complex, 2> B_cpy = B_gpu;
	BOOST_REQUIRE_CLOSE( B_cpy[1][0].real() , -0.72562939983295538 , 0.001);
	BOOST_REQUIRE_CLOSE( B_cpy[1][0].imag() ,  0.046772461520104877, 0.001);
}

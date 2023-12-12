// Copyright 2019-2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS dot"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/blas/dot.hpp>
#include <multi/array.hpp>

#include<complex>
#include<numeric>
#include<type_traits>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(blas_dot_context) {
	multi::array<float, 1> const x = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	blas::context const ctxt;
	{
		auto  res = +blas::dot(&ctxt, x, y);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.F) );
	}
	{
		float res = +blas::dot(&ctxt, x, y);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.F) );
	}
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context) {
	multi::array<float, 1> const x = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	auto res = +blas::dot(x, y);
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param) {
	multi::array<float, 1> const x = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	float res = NAN;
	blas::dot(x, y, multi::array_ref<float, 0>(res));
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.0F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex) {  // if you get a segfaut here, your system may require -DRETURN_BY_STACK
	using complex = std::complex<double>;
	multi::array<complex, 1> const x = { {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0} };  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = { {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0} };  // NOLINT(readability-identifier-length) BLAS naming
	complex res{0.0, 0.0};
	blas::dot(x, y, res);
	// an isolated error here might mean that the dot and nrm2 interface for the BLAS library is not detected properly
	BOOST_REQUIRE_EQUAL( real(res) , real(std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return alpha*std::conj(omega);})) );
	BOOST_REQUIRE_EQUAL( imag(res) , imag(std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return alpha*std::conj(omega);})) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C) {
	using complex = std::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
	complex res{0.0, 0.0};
	blas::dot(blas::C(x), y, res);
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
}

#if defined(CUDA_FOUND) and CUDA_FOUND
#include<thrust/complex.h>
BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C_thrust) {
	using complex = thrust::complex<double>; complex const I{0.0, 1.0};
	multi::array<complex, 1> const A = {1.0 + 0.*I, 2.0 + 0.0*I, 3.0 + 0.0*I};
	multi::array<complex, 1> const B = {1.0 + 0.*I, 2.0 + 2.0*I, 3.0 + 0.0*I};

	complex C;
	blas::dot(blas::C(A), B, C);
	BOOST_REQUIRE( C == std::inner_product(begin(A), end(A), begin(B), complex{0.0, 0.0}, std::plus<>{}, [](auto& a, auto& b){return conj(a)*b;}) );
}
#endif

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided) {
	multi::array<double, 2> const CA = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0}
	};
	double res = std::numeric_limits<double>::quiet_NaN();
	blas::dot_n(begin(CA[1]), size(CA[1]), begin(CA[2]), &res);
	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.0) );

	double const res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided_context) {
	multi::array<double, 2> const CA = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0}
	};
	double res = std::numeric_limits<double>::quiet_NaN();
	blas::context ctxt;
	blas::dot_n(&ctxt, begin(CA[1]), size(CA[1]), begin(CA[2]), &res);
	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.0) );

	double const res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_1d_real) {
	multi::array<float, 1> const x = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming

	using blas::dot;
	BOOST_REQUIRE( 14.0 == dot(x, y) );
	BOOST_REQUIRE( dot(x, y) == 14.0F );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_real) {
	multi::array<double, 2> const cA = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0}
	};
	{
		double const res = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( res == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0) );
	}
	{
		double res = NAN;
		blas::dot(cA[1], cA[2], res);
		BOOST_REQUIRE( res == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0) );
	}
	{
		double res = NAN;
		auto res2 = blas::dot(cA[1], cA[2], res);
		BOOST_REQUIRE( res == res2 );
	}
	 {
		double const res = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( res == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0) );
		BOOST_REQUIRE( blas::dot(cA[1], cA[2]) == blas::dot(cA[2], cA[1]) );
	}
}

BOOST_AUTO_TEST_CASE(inq_case) {
	multi::array<double, 1> x(multi::extensions_t<1>{multi::iextension{10}}, +1.0);  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> const y(multi::extensions_t<1>{multi::iextension{10}}, -1.0);  // NOLINT(readability-identifier-length) BLAS naming

	using blas::dot;
	using blas::hermitized;
	using blas::conj;

	auto res = dot(x, y);
	auto res2 = dot(hermitized(x), y);

	BOOST_REQUIRE(res == res2);

	auto res3 = dot(blas::conj(x), y);  // conjugation doesn't do anything for real array
	BOOST_REQUIRE(res3 == res);

	auto d_arr = dot(blas::C(x), y);
	BOOST_REQUIRE(d_arr == res);

	static_assert( not std::is_same<decltype(d_arr), double>{} );

	using blas::C;
	double const d_doub = dot(C(x), y);

	BOOST_REQUIRE( d_doub == d_arr );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex) {
	namespace blas = multi::blas;

	using complex = std::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 +     I,  2.0 + 3.0*I,  3.0 + 2.0*I,  4.0 - 9.0*I},
		{5.0 + 2.0*I,  6.0 + 6.0*I,  7.0 + 2.0*I,  8.0 - 3.0*I},
		{9.0 + 1.0*I, 10.0 + 9.0*I, 11.0 + 1.0*I, 12.0 + 2.0*I}
	};
	{
		complex c{0.0, 0.0}; blas::dot(A[1], A[2], c);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}) );
	}
	{
		complex c = blas::dot(A[1], A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}) );
	}
	{
		complex c = blas::dot(A[1], blas::C(A[2]));  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) {return alpha*conj(omega);}) );
	}
	{
		complex c = blas::dot(blas::C(A[1]), A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) {return conj(alpha)*omega;}) );
	}
	{
		complex c = blas::dot(blas::conj(A[1]), A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) {return conj(alpha)*omega;}) );
	}
	{
		complex c = blas::dot(blas::C(A[1]), A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) {return conj(alpha)*omega;}) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_conj_second) {
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	// using T = complex;
	using Alloc =  std::allocator<complex>;  // thrust::cuda::allocator<complex>;
	complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length)

	multi::array<complex, 2, Alloc> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {  // NOLINT(readability-identifier-length) BLAS naming
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({2, 2}, {3.0, 0.0});  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto CC = C;
		{
			auto const [is, js] = C.extensions();
			for(auto i : is) {
				for(auto j : js) {
					C[i][j] *= 0.0;
					for(auto k : B.extension()) {C[i][j] += A[i][k]*conj(B[k][j]);}  // NOLINT(altera-unroll-loops)
				}
			}
		}
		// TODO(correaa) MKL gives an error here
		// unknown location(0): fatal error: in "cublas_one_gemv_complex_conjtrans_zero": memory access violation at address: 0x00000007: no mapping at fault address
		{
			std::transform(begin(A), end(A), begin(CC), begin(CC), [BT = transposed(B)](auto const& Ar, auto&& Cr) {
				return std::transform(
					begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& Bc, auto&& Ce) {
						return 1.0*blas::dot(Ar, blas::C(Bc)) + 0.0*Ce;  // NOLINT(fuchsia-default-arguments-calls)
					}
				), std::forward<decltype(Cr)>(Cr);
			});
		}
		BOOST_TEST_REQUIRE( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST_REQUIRE( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST_REQUIRE( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST_REQUIRE( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
}

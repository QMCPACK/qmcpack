// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/test/unit_test.hpp>

#include <boost/multi/adaptors/blas/dot.hpp>
#include <boost/multi/array.hpp>

#include <complex>
#include <numeric>
#include <type_traits>

namespace multi = boost::multi;
namespace blas  = multi::blas;

BOOST_AUTO_TEST_CASE(blas_dot_context_double) {
	multi::array<double, 1> const x = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> const y = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming

	blas::context const ctxt;

	auto res1 = +blas::dot(&ctxt, x, y);
	BOOST_TEST( res1 == std::inner_product(begin(x), end(x), begin(y), 0.0) );

	auto const res2 = +blas::dot(&ctxt, x, y);
	BOOST_TEST( res2 == std::inner_product(begin(x), end(x), begin(y), 0.0) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_double) {
	multi::array<double, 1> const x = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> const y = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming

	auto res = +blas::dot(x, y);

	BOOST_TEST( res == std::inner_product(begin(x), end(x), begin(y), 0.0) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_double) {
	multi::array<double, 1> const x = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> const y = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming

	double res = NAN;

	blas::dot(x, y, multi::array_ref<double, 0>(res));
	BOOST_TEST( res == std::inner_product(begin(x), end(x), begin(y), 0.0) );
}

// float uses of dot are disabled because of a bug in Apple Accelerate BLAS, https://fortran-lang.discourse.group/t/how-many-blas-libraries-have-this-error/4454/23, https://forums.developer.apple.com/forums/thread/717757
BOOST_AUTO_TEST_CASE(blas_dot_context_float, *boost::unit_test::disabled()) {
	multi::array<float, 1> const x = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming

	blas::context const ctxt;

	auto res1 = +blas::dot(&ctxt, x, y);
	BOOST_TEST( res1 == std::inner_product(begin(x), end(x), begin(y), 0.0F) );

	auto const res2 = +blas::dot(&ctxt, x, y);
	BOOST_TEST( res2 == std::inner_product(begin(x), end(x), begin(y), 0.0F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_float, *boost::unit_test::disabled()) {
	multi::array<float, 1> const x = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming

	auto res = +blas::dot(x, y);

	BOOST_TEST( res == std::inner_product(begin(x), end(x), begin(y), 0.0F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_float, *boost::unit_test::disabled()) {
	multi::array<float, 1> const x   = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y   = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	float                        res = NAN;
	blas::dot(x, y, multi::array_ref<float, 0>(res));
	BOOST_TEST( res == std::inner_product(begin(x), end(x), begin(y), 0.0F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_double) {  // if you get a segfaut here, your system may require -DRETURN_BY_STACK
	using complex = std::complex<double>;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const x = {
		{1.0, 0.0},
		{2.0, 0.0},
		{3.0, 0.0},
	};
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = {
		{1.0, 0.0},
		{2.0, 0.0},
		{3.0, 0.0},
	};  // NOLINT(readability-identifier-length) BLAS naming
	complex res{0.0, 0.0};
	blas::dot(x, y, res);
	// an isolated error here might mean that the dot and nrm2 interface for the BLAS library is not detected properly
	BOOST_REQUIRE_EQUAL(real(res), real(std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) { return alpha * std::conj(omega); })));
	BOOST_REQUIRE_EQUAL(imag(res), imag(std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) { return alpha * std::conj(omega); })));
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_float) {  // if you get a segfaut here, your system may require -DRETURN_BY_STACK
	using complex = std::complex<float>;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const x = {
		{1.0F, 0.0F},
		{2.0F, 0.0F},
		{3.0F, 0.0F},
	};
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = {
		{1.0F, 0.0F},
		{2.0F, 0.0F},
		{3.0F, 0.0F},
	};  // NOLINT(readability-identifier-length) BLAS naming
	complex res{0.0F, 0.0F};
	blas::dot(x, y, res);

	// // an isolated error here might mean that the dot and nrm2 interface for the BLAS library is not detected properly
	BOOST_REQUIRE_EQUAL(real(res), real(std::inner_product(begin(x), end(x), begin(y), complex{0.0F, 0.0F}, std::plus<>{}, [](auto const& alpha, auto const& omega) { return alpha * std::conj(omega); })));
	BOOST_REQUIRE_EQUAL(imag(res), imag(std::inner_product(begin(x), end(x), begin(y), complex{0.0F, 0.0F}, std::plus<>{}, [](auto const& alpha, auto const& omega) { return alpha * std::conj(omega); })));
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C) {
	using complex = std::complex<double>;
	auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> const x = {1.0 + 0.0 * I, 2.0 + 0.0 * I, 3.0 + 0.0 * I};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = {1.0 + 0.0 * I, 2.0 + 2.0 * I, 3.0 + 0.0 * I};  // NOLINT(readability-identifier-length) BLAS naming

	complex res{0.0, 0.0};
	blas::dot(blas::C(x), y, res);
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) { return conj(alpha) * omega;}) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C_float) {
	using complex = std::complex<float>;
	auto const I  = complex{0.0F, 1.0F};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> const x = {1.0F + 0.0F * I, 2.0F + 0.0F * I, 3.0F + 0.0F * I};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = {1.0F + 0.0F * I, 2.0F + 2.0F * I, 3.0F + 0.0F * I};  // NOLINT(readability-identifier-length) BLAS naming

	complex res{0.0F, 0.0F};
	blas::dot(blas::C(x), y, res);
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0F, 0.0F}, std::plus<>{}, [](auto const& alpha, auto const& omega) { return conj(alpha) * omega;}) );
}

#if defined(CUDA_FOUND) and CUDA_FOUND
#include <thrust/complex.h>
BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C_thrust) {
	using complex = thrust::complex<double>;
	auto const I  = complex{0.0, 1.0};

	multi::array<complex, 1> const A = {1.0 + 0. * I, 2.0 + 0.0 * I, 3.0 + 0.0 * I};
	multi::array<complex, 1> const B = {1.0 + 0. * I, 2.0 + 2.0 * I, 3.0 + 0.0 * I};

	complex C;
	blas::dot(blas::C(A), B, C);
	BOOST_REQUIRE( C == std::inner_product(begin(A), end(A), begin(B), complex{0.0, 0.0}, std::plus<>{}, [](auto& a, auto& b){ return conj(a) * b;}) );
}
#endif

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided_double) {
	multi::array<double, 2> const CA = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	double res = std::numeric_limits<double>::quiet_NaN();
	blas::dot_n(begin(CA[1]), size(CA[1]), begin(CA[2]), &res);
	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.0) );

	double const res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided_float, *boost::unit_test::disabled()) {
	multi::array<float, 2> const CA = {
		{1.0F,  2.0F,  3.0F,  4.0F},
		{5.0F,  6.0F,  7.0F,  8.0F},
		{9.0F, 10.0F, 11.0F, 12.0F},
	};
	auto res = std::numeric_limits<float>::quiet_NaN();
	blas::dot_n(begin(CA[1]), size(CA[1]), begin(CA[2]), &res);
	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.0F) );

	double const res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided_context) {
	multi::array<double, 2> const CA = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};
	double        res = std::numeric_limits<double>::quiet_NaN();
	blas::context ctxt;
	blas::dot_n(&ctxt, begin(CA[1]), size(CA[1]), begin(CA[2]), &res);
	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.0) );

	double const res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided_context_float, *boost::unit_test::disabled()) {
	multi::array<float, 2> const CA = {
		{1.0F,  2.0F,  3.0F,  4.0F},
		{5.0F,  6.0F,  7.0F,  8.0F},
		{9.0F, 10.0F, 11.0F, 12.0F},
	};
	float res = std::numeric_limits<double>::quiet_NaN();

	blas::context ctxt;
	blas::dot_n(&ctxt, begin(CA[1]), size(CA[1]), begin(CA[2]), &res);

	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.0) );

	float const res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_1d_real_double) {
	multi::array<double, 1> const x = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> const y = {1.0, 2.0, 3.0};  // NOLINT(readability-identifier-length) BLAS naming

	using blas::dot;
	BOOST_TEST( 14.0 == dot(x, y) );
	BOOST_TEST( dot(x, y) == 14.0F );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_1d_real_float, *boost::unit_test::disabled()) {
	multi::array<float, 1> const x = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1.0F, 2.0F, 3.0F};  // NOLINT(readability-identifier-length) BLAS naming

	using blas::dot;
	BOOST_TEST( 14.0F == dot(x, y) );
	BOOST_TEST( dot(x, y) == 14.0F );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_real_double) {
	multi::array<double, 2> const cA = {
		{1.0,  2.0,  3.0,  4.0},
		{5.0,  6.0,  7.0,  8.0},
		{9.0, 10.0, 11.0, 12.0},
	};

	double const res1 = blas::dot(cA[1], cA[2]);
	BOOST_REQUIRE( res1 == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0) );

	double res2 = NAN;
	blas::dot(cA[1], cA[2], res2);
	BOOST_REQUIRE( res2 == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0) );

	double       res_nan = NAN;
	double const res3    = blas::dot(cA[1], cA[2], res_nan);
	BOOST_REQUIRE( res3 == res2 );

	double const res4 = blas::dot(cA[1], cA[2]);
	BOOST_REQUIRE( res4 == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0) );
	BOOST_REQUIRE( blas::dot(cA[1], cA[2]) == blas::dot(cA[2], cA[1]) );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_real_float, *boost::unit_test::disabled()) {
	multi::array<float, 2> const cA = {
		{1.0F,  2.0F,  3.0F,  4.0F},
		{5.0F,  6.0F,  7.0F,  8.0F},
		{9.0F, 10.0F, 11.0F, 12.0F},
	};

	float const res1 = blas::dot(cA[1], cA[2]);
	BOOST_REQUIRE( res1 == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0F) );

	float res2 = NAN;
	blas::dot(cA[1], cA[2], res2);
	BOOST_REQUIRE( res2 == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0F) );

	float       res_nan = NAN;
	float const res3    = blas::dot(cA[1], cA[2], res_nan);

	BOOST_REQUIRE( res3 == res2 );

	float const res4 = blas::dot(cA[1], cA[2]);
	BOOST_REQUIRE( res4 == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.0F) );
	BOOST_REQUIRE( blas::dot(cA[1], cA[2]) == blas::dot(cA[2], cA[1]) );
}

BOOST_AUTO_TEST_CASE(inq_case) {
	multi::array<double, 1> const x(multi::extensions_t<1>{multi::iextension{10}}, +1.0);  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> const y(multi::extensions_t<1>{multi::iextension{10}}, -1.0);  // NOLINT(readability-identifier-length) BLAS naming

	using blas::conj;
	using blas::dot;
	using blas::hermitized;

	auto res  = dot(x, y);
	auto res2 = dot(hermitized(x), y);

	BOOST_REQUIRE(res == res2);

	auto res3 = dot(blas::conj(x), y);  // conjugation doesn't do anything for real array
	BOOST_REQUIRE(res3 == res);

	auto d_arr = dot(blas::C(x), y);
	BOOST_REQUIRE(d_arr == res);

	static_assert(!std::is_same<decltype(d_arr), double>{});

	using blas::C;
	double const d_doub = dot(C(x), y);

	BOOST_REQUIRE( d_doub == d_arr );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex_double) {
	namespace blas = multi::blas;

	using complex = std::complex<double>;

	complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const A = {
		{      1.0 + I,  2.0 + 3.0 * I,  3.0 + 2.0 * I,  4.0 - 9.0 * I},
		{5.0 + 2.0 * I,  6.0 + 6.0 * I,  7.0 + 2.0 * I,  8.0 - 3.0 * I},
		{9.0 + 1.0 * I, 10.0 + 9.0 * I, 11.0 + 1.0 * I, 12.0 + 2.0 * I},
	};

	auto c1 = complex{0.0, 0.0};
	blas::dot(A[1], A[2], c1);
	BOOST_TEST_REQUIRE( c1 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}) );

	auto const c2 = +blas::dot(A[1], A[2]);
	BOOST_TEST_REQUIRE( c2 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}) );

	complex const c3 = blas::dot(A[1], A[2]);
	BOOST_TEST_REQUIRE( c3 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}) );

	complex const c4 = blas::dot(A[1], blas::C(A[2]));
	BOOST_TEST_REQUIRE( c4 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) { return alpha * conj(omega);}) );

	complex const c5 = blas::dot(blas::C(A[1]), A[2]);
	BOOST_TEST_REQUIRE( c5 == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) { return conj(alpha) * omega;}) );

	complex const c6 = blas::dot(blas::conj(A[1]), A[2]);
	BOOST_TEST_REQUIRE( c6 == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) { return conj(alpha) * omega;}) );

	complex const c7 = blas::dot(blas::C(A[1]), A[2]);
	BOOST_TEST_REQUIRE( c7 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0, 0.0}, std::plus<>{}, [](auto alpha, auto omega) { return conj(alpha) * omega;}) );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex_float) {
	namespace blas = multi::blas;

	using complex = std::complex<float>;

	complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2> const A = {
		{1.0F + 1.0F * I,  2.0F + 3.0F * I,  3.0F + 2.0F * I,  4.0F - 9.0F * I},
		{5.0F + 2.0F * I,  6.0F + 6.0F * I,  7.0F + 2.0F * I,  8.0F - 3.0F * I},
		{9.0F + 1.0F * I, 10.0F + 9.0F * I, 11.0F + 1.0F * I, 12.0F + 2.0F * I},
	};

	auto c1 = complex{0.0F, 0.0F};
	blas::dot(A[1], A[2], c1);
	BOOST_TEST_REQUIRE( c1 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0F, 0.0F}) );

	auto const c2 = +blas::dot(A[1], A[2]);
	BOOST_TEST_REQUIRE( c2 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0F, 0.0F}) );

	complex const c3 = blas::dot(A[1], A[2]);
	BOOST_TEST_REQUIRE( c3 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0F, 0.0F}) );

	complex const c4 = blas::dot(A[1], blas::C(A[2]));
	BOOST_TEST_REQUIRE( c4 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0F, 0.0F}, std::plus<>{}, [](auto alpha, auto omega) { return alpha * conj(omega);}) );

	complex const c5 = blas::dot(blas::C(A[1]), A[2]);
	BOOST_TEST_REQUIRE( c5 == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0F, 0.0F}, std::plus<>{}, [](auto alpha, auto omega) { return conj(alpha) * omega;}) );

	complex const c6 = blas::dot(blas::conj(A[1]), A[2]);
	BOOST_TEST_REQUIRE( c6 == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0F, 0.0F}, std::plus<>{}, [](auto alpha, auto omega) { return conj(alpha) * omega;}) );

	complex const c7 = blas::dot(blas::C(A[1]), A[2]);
	BOOST_TEST_REQUIRE( c7 == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.0F, 0.0F}, std::plus<>{}, [](auto alpha, auto omega) { return conj(alpha) * omega;}) );
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_conj_second_double) {
	namespace blas = multi::blas;

	using complex = std::complex<double>;
	using Alloc   = std::allocator<complex>;  // thrust::cuda::allocator<complex>;

	auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};

	multi::array<complex, 2, Alloc> C({2, 2}, {3.0, 0.0});  // NOLINT(readability-identifier-length) conventional BLAS naming

	auto CC = C;

	auto const [is, js] = C.extensions();
	std::for_each(is.begin(), is.end(), [&, js = js](auto ii) {
		std::for_each(js.begin(), js.end(), [&](auto jj) {
			C[ii][jj] *= 0.0;
			std::for_each(B.extension().begin(), B.extension().end(), [&](auto kk) {
				C[ii][jj] += A[ii][kk] * conj(B[kk][jj]);
			});
		});
	});

	// TODO(correaa) MKL gives an error here
	// unknown location(0): fatal error: in "cublas_one_gemv_complex_conjtrans_zero": memory access violation at address: 0x00000007: no mapping at fault address

	std::transform(begin(A), end(A), begin(CC), begin(CC), [BT = transposed(B)](auto const& Ar, auto&& Cr) {
		return std::transform(
			       begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& Bc, auto const& Ce) {
				       return std::complex<double>{1.0, 0.0} * blas::dot(Ar, blas::C(Bc)) + 0.0 * Ce;
			       }
		       ),
		       std::forward<decltype(Cr)>(Cr);
	});

	BOOST_TEST_REQUIRE( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	BOOST_TEST_REQUIRE( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	BOOST_TEST_REQUIRE( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	BOOST_TEST_REQUIRE( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_conj_second_float) {
	namespace blas = multi::blas;

	using complex = std::complex<float>;
	using Alloc   = std::allocator<complex>;  // thrust::cuda::allocator<complex>;

	auto const I = complex{0.0F, 1.0F};  // NOLINT(readability-identifier-length)

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0F - 2.0F * I, 9.0F - 1.0F * I},
		{2.0F + 3.0F * I, 1.0F - 2.0F * I},
	};

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const B = {
		{3.0F - 4.0F * I, 19.0F - 1.0F * I},
		{1.0F + 5.0F * I,  8.0F - 8.0F * I},
	};

	multi::array<complex, 2, Alloc> C({2, 2}, {3.0F, 0.0F});  // NOLINT(readability-identifier-length) conventional BLAS naming

	auto CC = C;

	auto const [is, js] = C.extensions();
	std::for_each(is.begin(), is.end(), [&, js = js](auto ii) {
		std::for_each(js.begin(), js.end(), [&](auto jj) {
			C[ii][jj] *= 0.0F;
			std::for_each(B.extension().begin(), B.extension().end(), [&](auto kk) {
				C[ii][jj] += A[ii][kk] * conj(B[kk][jj]);
			});
		});
	});

	// TODO(correaa) MKL gives an error here
	// unknown location(0): fatal error: in "cublas_one_gemv_complex_conjtrans_zero": memory access violation at address: 0x00000007: no mapping at fault address

	std::transform(begin(A), end(A), begin(CC), begin(CC), [BT = transposed(B)](auto const& Ar, auto&& Cr) {
		return std::transform(
			       begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& Bc, auto const& Ce) {
				       return complex{1.0F, 0.0F} * blas::dot(Ar, blas::C(Bc)) + 0.0F * Ce;
			       }
		       ),
		       std::forward<decltype(Cr)>(Cr);
	});

	BOOST_TEST_REQUIRE( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	BOOST_TEST_REQUIRE( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	BOOST_TEST_REQUIRE( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	BOOST_TEST_REQUIRE( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
}

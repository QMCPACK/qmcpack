// Copyright 2023-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/cuda/cublas.hpp>

#include <boost/multi/adaptors/blas/axpy.hpp>
#include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/adaptors/blas/gemv.hpp>
#include <boost/multi/adaptors/blas/nrm2.hpp>

#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/complex.h>

#include <numeric>  // for std::iota

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>

template<class T>
void test_gemv_complex(double tol) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<T>;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<complex, 2> const M_gpu = {
		{ {T( 9.0), T(0.0)}, {T(24.0), T(0.0)}, {T(30.0), T(0.0)}, {T(9.0), T(0.0)}},
		{ {T( 4.0), T(0.0)}, {T(10.0), T(0.0)}, {T(12.0), T(0.0)}, {T(7.0), T(0.0)}},
		{{T(14.0), T(0.0)}, {T(16.0), T(0.0)}, {T(36.0), T(0.0)}, {T(1.0), T(0.0)}},
	};

	multi::thrust::cuda::array<complex, 1> const X_gpu = {
		{T(1.1), T(0.0)},
		{T(2.1), T(0.0)},
		{T(3.1), T(0.0)},
		{T(4.1), T(0.0)}
	};

	multi::thrust::cuda::array<complex, 1> Y_gpu = {
		{T(4.0), T(0.0)},
		{T(5.0), T(0.0)},
		{T(6.0), T(0.0)}
	};

	blas::gemv(/*alpha*/ T(1.1), M_gpu, X_gpu, /*beta*/ T(1.2), Y_gpu);  // y = a*M*x + b*y

	multi::array<complex, 1> const Y_copy = Y_gpu;

	using blas::operators::operator-;
	BOOST_TEST(+blas::nrm2(Y_copy -
		multi::array<complex, 1>{
			{T(214.02), T(0.0)},
			{T(106.43), T(0.0)},
			{T(188.37), T(0.0)}
		}
	) < tol);
}

template<class T>
void test_gemv_real(double tol) {
	namespace blas = multi::blas;
	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<T, 2> const M_gpu = {
		{T( 9.0), T(24.0), T(30.0), T(9.0)},
		{T( 4.0), T(10.0), T(12.0), T(7.0)},
		{T(14.0), T(16.0), T(36.0), T(1.0)},
	};

	multi::thrust::cuda::array<T, 1> const X_gpu = {T(1.1), T(2.1), T(3.1), T(4.1)};

	multi::thrust::cuda::array<T, 1> Y_gpu = {T(4.0), T(5.0), T(6.0)};

	blas::gemv(/*alpha*/ T(1.1), M_gpu, X_gpu, /*beta*/ T(1.2), Y_gpu);  // y = a*M*x + b*y

	multi::array<T, 1> const Y_copy = Y_gpu;

	using blas::operators::operator-;
	BOOST_TEST(+blas::nrm2(Y_copy - multi::array<T, 1>{T(214.02), T(106.43), T(188.37)}) < tol);
}

template<class T>
void test_gemm_real() {
	namespace blas = multi::blas;
	{
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::thrust::cuda::array<T, 2> const a = {
			{T(1.0), T(9.0)},
			{T(2.0), T(1.0)},
		};

		multi::thrust::cuda::array<T, 2> c({2, 2}, T(9999.0));  // NOLINT(readability-identifier-length) conventional BLAS naming
		blas::gemm(T(1.0), a, a, T(0.0), c);  // c=aa

		multi::array<T, 2> const c_copy = c;
		BOOST_TEST( c_copy[0][0] == T(19.0) );
		BOOST_TEST( c_copy[0][1] == T(18.0) );
		BOOST_TEST( c_copy[1][0] == T( 4.0) );
		BOOST_TEST( c_copy[1][1] == T(19.0) );
	}
	{
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::thrust::cuda::array<T, 2> const a = {
			{T(1.0), T(9.0)},
			{T(2.0), T(1.0)},
		};

		auto const c = +blas::gemm(T(1.0), a, a);  // c=aa

		multi::array<T, 2> const c_copy = c;
		BOOST_TEST( c_copy[0][0] == T(19.0) );
		BOOST_TEST( c_copy[0][1] == T(18.0) );
		BOOST_TEST( c_copy[1][0] == T( 4.0) );
		BOOST_TEST( c_copy[1][1] == T(19.0) );
	}
}

template<class T>
void test_gemm_complex() {
	namespace blas = multi::blas;

	using complex = thrust::complex<T>;
	complex const I{T(0.0), T(1.0)};  // NOLINT(readability-identifier-length) imaginary unit

	{
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::thrust::cuda::array<complex, 2> const a = {
			{T(1.0) - T(2.0) * I, T(9.0) - T(1.0) * I},
			{T(2.0) + T(3.0) * I, T(1.0) - T(2.0) * I},
		};

		multi::thrust::cuda::array<complex, 2> c({2, 2}, {T(9999.0), T(0.0)});  // NOLINT(readability-identifier-length) conventional BLAS naming
		blas::gemm({T(1.0), T(0.0)}, a, a, {T(0.0), T(0.0)}, c);  // c=aa†, c†=aa†

		multi::array<complex, 2> const c_copy = c;
		BOOST_TEST( c_copy[1][0] == T(16.0) - T( 2.0)*I );
		BOOST_TEST( c_copy[0][1] == T(14.0) - T(38.0)*I );
	}
	{
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::thrust::cuda::array<complex, 2> const a = {
			{T(1.0) - T(2.0) * I, T(9.0) - T(1.0) * I},
			{T(2.0) + T(3.0) * I, T(1.0) - T(2.0) * I},
		};

		auto const c = +blas::gemm(complex{T(1.0), T(0.0)}, a, a);  // c=aa†, c†=aa†

		multi::array<complex, 2> const c_copy = c;
		BOOST_TEST( c_copy[1][0] == T(16.0) - T( 2.0)*I );
		BOOST_TEST( c_copy[0][1] == T(14.0) - T(38.0)*I );
	}
	{
		// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
		multi::thrust::cuda::array<complex, 2> const a = {
			{T(1.0) - T(2.0) * I, T(9.0) - T(1.0) * I},
			{T(2.0) + T(3.0) * I, T(1.0) - T(2.0) * I},
		};

		multi::thrust::cuda::array<complex, 2> c({2, 2}, {T(0.0), T(0.0)});
		c += blas::gemm(complex{T(1.0), T(0.0)}, a, a);  // c=aa†, c†=aa†

		multi::array<complex, 2> const c_copy = c;
		BOOST_TEST( c_copy[1][0] == T(16.0) - T( 2.0)*I );
		BOOST_TEST( c_copy[0][1] == T(14.0) - T(38.0)*I );
	}
}

int main() {
	test_gemv_real<float> (1e-4);
	test_gemv_real<double>(1e-13);

	test_gemv_complex<float> (1e-4);
	test_gemv_complex<double>(1e-13);

	test_gemm_real<float>();
	test_gemm_real<double>();

	test_gemm_complex<float>();
	test_gemm_complex<double>();

	// Chris
	{
		int m   = 12;
		int n   = 10;
		int ell = 8;

		multi::array<double, 3> a({m, n, ell});
		multi::array<double, 2> b({n, ell});

		std::iota(a.elements().begin(), a.elements().end(), 20.0);
		std::iota(b.elements().begin(), b.elements().end(), 30.0);

		multi::array<double, 1> c_gold(m, 0.0);

		for(int k = 0; k != m; ++k) {
			for(int j = 0; j != n; ++j) {
				for(int i = 0; i != ell; ++i) {
					c_gold[k] += a[k][j][i] * b[j][i];
				}
			}
		}

		multi::array<double, 1> c_flat(m, 0.0);

		for(int k = 0; k != m; ++k) {
			for(int ji = 0; ji != a[k].elements().size(); ++ji) {
				c_flat[k] += a[k].elements()[ji] * b.elements()[ji];
			}
		}

		BOOST_TEST( c_gold == c_flat );

		// gpu::run(m, n, [...] (auto k, auto j) {
		//      for (int i =0; i <l ; ++i) {
		//        c[k] +=a[k][j][i]) * b[j][i];
		//        }
		//  });
	}

	return boost::report_errors();
}

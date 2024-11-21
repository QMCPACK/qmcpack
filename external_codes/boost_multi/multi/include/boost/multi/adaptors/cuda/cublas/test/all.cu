// Copyright 2023 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS all"
// #include <boost/test/unit_test.hpp>

#include <boost/multi/adaptors/cuda/cublas.hpp>

#include <boost/multi/adaptors/blas/asum.hpp>
#include <boost/multi/adaptors/blas/axpy.hpp>
#include <boost/multi/adaptors/blas/copy.hpp>
#include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/adaptors/blas/gemv.hpp>
#include <boost/multi/adaptors/blas/nrm2.hpp>
#include <boost/multi/adaptors/blas/scal.hpp>
#include <boost/multi/adaptors/blas/swap.hpp>
#include <boost/multi/adaptors/blas/trsm.hpp>

#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/complex.h>

#include <numeric>
#include <thrust/inner_product.h>
#include <thrust/transform_reduce.h>

namespace multi = boost::multi;

using complex = thrust::complex<double>;

template<class T = complex, class Alloc = std::allocator<T>>
auto generate_ABx() {
	complex const             I{ 0.0, 1.0 };
	multi::array<T, 1, Alloc> x = { 1.0 + I * 0.0, 2.0 + I * 0.0, 3.0 + I * 0.0, 4.0 + I * 0.0 };

	multi::array<complex, 2, Alloc> A = {
		{ 1.0 + I * 0.0,  2.0 + I * 0.0,  3.0 + I * 0.0,  4.0 + I * 0.0},
		{ 5.0 + I * 0.0,  6.0 + I * 0.0,  7.0 + I * 0.0,  8.0 + I * 0.0},
		{ 9.0 + I * 0.0, 10.0 + I * 0.0, 11.0 + I * 0.0, 12.0 + I * 0.0},
		{13.0 + I * 0.0, 14.0 + I * 0.0, 15.0 + I * 0.0, 16.0 + I * 0.0},
	};

	multi::array<complex, 2, Alloc> B = {
		{ 1.0 + I * 0.0,  2.0 + I * 0.0,  3.0 + I * 0.0,  4.0 + I * 0.0},
		{ 5.0 + I * 0.0,  6.0 + I * 0.0,  7.0 + I * 0.0,  8.0 + I * 0.0},
		{ 9.0 + I * 0.0, 10.0 + I * 0.0, 11.0 + I * 0.0, 12.0 + I * 0.0},
		{13.0 + I * 0.0, 14.0 + I * 0.0, 15.0 + I * 0.0, 16.0 + I * 0.0},
	};

	return std::make_tuple(std::move(x), std::move(A), std::move(B));
}

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

#define BOOST_REQUIRE_CLOSE(X, Y, ToL) BOOST_TEST( std::abs( (X) - (Y) ) < (ToL) )

int main() {
BOOST_AUTO_TEST_CASE(cublas_scal_complex_column) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	{
		using T        = complex;
		auto [x, A, B] = generate_ABx<T, thrust::cuda::allocator<T>>();
		auto const s   = 2.0 + I * 3.0;
		blas::scal(s, x);  // x_i <- s*x_i

		{
			auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex>>();
			auto xx           = +x2;
			blas::scal(s, xx);
			BOOST_TEST(xx == x);
		}
		{
			auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex>>();
			using blas::operators::operator*=;
			x2 *= s;
			BOOST_TEST(x == x2);
		}
		{
			auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex>>();
			thrust::transform(x2.begin(), x2.end(), x2.begin(), [s] __device__(T & e) { return s * e; });

			BOOST_TEST(x == x2);
		}
		{
			auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex>>();
			thrust::for_each(x2.begin(), x2.end(), [s] __device__(T & e) { return e *= s; });

			BOOST_TEST(x == x2);
		}
	}
}

BOOST_AUTO_TEST_CASE(cublas_copy_complex) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<T>;

	multi::array<T, 1, Alloc> const x = { 1.0 + I * 8.0, 2.0 + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };
	multi::array<T, 1, Alloc>       y = { 1.0 + I * 9.0, 2.0 + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };

	blas::copy(x, y);
	BOOST_TEST( static_cast<complex>(y[0]) == 1.0 + I*8.0 );
	{
		thrust::copy(begin(x), end(x), begin(y));
		BOOST_TEST( static_cast<complex>(y[0]) == 1.0 + I*8.0 );
	}
	{
		blas::copy_n(x.begin(), x.size(), y.begin());
		BOOST_TEST( static_cast<complex>(y[0]) == 1.0 + I*8.0 );
	}
	{
		y() = blas::copy(x);
		BOOST_TEST( static_cast<complex>(y[0]) == 1.0 + I*8.0 );
	}
	{
		multi::array<T, 1, Alloc> yy = blas::copy(x);
		BOOST_TEST( static_cast<complex>(yy[0]) == 1.0 + I*8.0 );
	}
	{
		y = blas::copy(x);
		BOOST_TEST( static_cast<complex>(y[0]) == 1.0 + I*8.0 );
	}
	{
		{
			using blas::operators::operator<<;
			y << x;
			//  BOOST_TEST(( static_cast<complex>(y[0]) == 1.0 + I*8.0 ));  // this can't be used with a free operator<<
		}
		BOOST_TEST(( static_cast<complex>(y[0]) == 1.0 + I*8.0 ));  // this can't be used with a free operator<<
	}
}

#if 1
BOOST_AUTO_TEST_CASE(cublas_swap_complex) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<T>;

	multi::array<T, 1, Alloc> x = { 1.0 + I * 8.0, 2.0 + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };
	multi::array<T, 1, Alloc> y = { 1.0 + I * 9.0, 2.0 + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };

	blas::swap(x, y);
	BOOST_TEST( static_cast<complex>(x[0]) == 1.0 + I*9.0 );
	{
		thrust::swap_ranges(begin(x), end(x), begin(y));
		thrust::swap_ranges(begin(x), end(x), begin(y));
		BOOST_TEST( static_cast<complex>(x[0]) == 1.0 + I*9.0 );
	}
	{
		using blas::operator^;
		(x ^ y);
		(x ^ y);
		BOOST_TEST( static_cast<complex>(x[0]) == 1.0 + I*9.0 );
	}
}

BOOST_AUTO_TEST_CASE(cublas_asum_complex_column) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<T, 1, Alloc> const x = { 1.0 + I * 8.0, 2.0 + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };

	double res;
	blas::asum_n(x.begin(), x.size(), &res);
	{
		double res2;
		res2 = blas::asum(x);
		BOOST_TEST( res == res2 );
	}
	{
		double res2 = blas::asum(x);
		BOOST_TEST( res == res2 );
	}
	{
		auto res2 = std::transform_reduce(
			x.begin(), x.end(), double{}, std::plus<>{}, [](T const& e) { return std::abs(e.real()) + std::abs(e.imag()); }
		);
		BOOST_TEST( res == res2 );
	}
	{
		auto res2 = thrust::transform_reduce(
			x.begin(), x.end(),
			[] __host__ __device__(T const& e) { return std::abs(e.real()) + std::abs(e.imag()); },
			double{}, thrust::plus<>{}
		);
		BOOST_TEST( res == res2 );
	}
	{
		multi::static_array<double, 0, thrust::cuda::allocator<double>> res2({}, 0.0);
		res2.assign(&blas::asum(x));
		res2 = blas::asum(x);
		BOOST_TEST(( res == static_cast<multi::static_array<double, 0, thrust::cuda::allocator<double>>::element_ref>(res2) ));
		BOOST_TEST(( res == static_cast<double>(res2) ));
		//  BOOST_TEST( res == res2 );
	}
	{
		multi::array<double, 0, thrust::cuda::allocator<double>> res2 = blas::asum(x);
		BOOST_TEST(( res == static_cast<multi::static_array<double, 0, thrust::cuda::allocator<double>>::element_ref>(res2) ));
		BOOST_TEST(( res == static_cast<double>(res2) ));
		//  BOOST_TEST( res == res2 );
	}
	{
		using blas::operators::operator==;
		using blas::operators::operator!=;
		BOOST_TEST( x != 0 );
		BOOST_TEST( not (x == 0) );
	}
	{
		using blas::operators::contains_nan;
		BOOST_TEST( not contains_nan(x) );
	}
	{
		using blas::operators::isfinite;
		using blas::operators::isinf;
		BOOST_TEST( isfinite(x) );
		BOOST_TEST( not isinf(x) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_asum_complex_nans) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<T, 1, Alloc> const x = { 1.0 + I * 8.0, std::numeric_limits<double>::quiet_NaN() + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };

	{
		using blas::operators::contains_nan;
		BOOST_TEST( contains_nan(x) );
	}
	{
		using blas::operators::operator==;
		using blas::operators::operator!=;
		BOOST_TEST( not (x != 0) );
		BOOST_TEST( not (x == 0) );
	}
	{
		using blas::operators::isfinite;
		using blas::operators::isinf;
		BOOST_TEST( not isfinite(x) );
		BOOST_TEST( not isinf(x) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_asum_complex_inf) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<T, 1, Alloc> const x = { 1.0 + I * 8.0, std::numeric_limits<double>::infinity() + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };

	// double res;
	{
		using blas::operators::contains_nan;
		BOOST_TEST( not contains_nan(x) );
	}
	{
		using blas::operators::operator==;
		using blas::operators::operator!=;
		BOOST_TEST(     (x != 0) );
		BOOST_TEST( not (x == 0) );
	}
	{
		using blas::operators::isfinite;
		using blas::operators::isinf;
		BOOST_TEST( not isfinite(x) );
		BOOST_TEST( isinf(x) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_nrm2_complex_column) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<T, 1, Alloc> const x = { 1.0 + I * 8.0, 2.0 + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };

	double res;
	blas::nrm2(x, res);
	{
		double res2;
		res2 = blas::nrm2(x);
		BOOST_TEST( res == res2 );
	}
	{
		auto res2 = +blas::nrm2(x);
		BOOST_TEST( res == res2 );
	}
	{
		auto res2 = sqrt(thrust::transform_reduce(
			x.begin(), x.end(),
			[] __host__ __device__(T const& e) { return thrust::norm(e); },
			double{}, thrust::plus<>{}
		));
		BOOST_TEST( res == res2 );
	}
	{
		multi::array<double, 0, thrust::cuda::allocator<double>> res2 = blas::nrm2(x);
		BOOST_TEST(( res == static_cast<double>(res2) ));
	}
}

BOOST_AUTO_TEST_CASE(cublas_dot_complex_column) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<T, 1, Alloc> const x = { 1.0 + I * 8.0, 2.0 + I * 6.0, 3.0 + I * 5.0, 4.0 + I * 3.0 };
	multi::array<T, 1, Alloc> const y = { 1.0 + I * 2.0, 2.0 + I * 3.0, 3.0 + I * 5.0, 4.0 + I * 7.0 };

	{
		T res;
		blas::dot(x, y, res);
		{
			complex res2;
			res2 = blas::dot(x, y);
			BOOST_TEST(res == res2);
		}
		{
			multi::array<complex, 0> res2(complex{ 1.0, 0.0 });
			res2 = blas::dot(x, y);
			BOOST_TEST( static_cast<complex>(res2) == res );
		}
		{
			using blas::operators::operator, ;
			auto res2 = +(x, y);
			BOOST_TEST(res == res2);
		}
		{
			auto res2 = +blas::dot(x, y);
			BOOST_TEST(res == res2);
		}
		{
			//  auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex> >();
			//  thrust::for_each(x2.begin(), x2.end(), [s] __device__ (T& e) {return e*=s;});
			auto res2 = thrust::inner_product(x.begin(), x.end(), y.begin(), T{});
			BOOST_TEST(res == res2);
		}
	}
	{
		T res;
		blas::dot(blas::C(x), y, res);
		{
			using blas::operators::operator, ;
			using blas::operators::operator*;
			auto res2 = +(*x, y);
			BOOST_TEST(res == res2);
		}
		{
			auto res2 = +blas::dot(blas::C(x), y);
			BOOST_TEST(res == res2);
		}
		{
			//  auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex> >();
			//  thrust::for_each(x2.begin(), x2.end(), [s] __device__ (T& e) {return e*=s;});
			auto res2 = thrust::inner_product(x.begin(), x.end(), y.begin(), T{}, thrust::plus<>{}, [] __device__(T const& t1, T const& t2) { return conj(t1) * t2; });
			BOOST_TEST(res == res2);
		}
	}
	{
		T res;
		blas::dot(x, blas::C(y), res);
		{
			using blas::operators::operator, ;
			auto res2 = +(x, blas::C(y));
			BOOST_TEST(res == res2);
		}
		{
			auto res2 = +blas::dot(x, blas::C(y));
			BOOST_TEST(res == res2);
		}
		{
			//  auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex> >();
			//  thrust::for_each(x2.begin(), x2.end(), [s] __device__ (T& e) {return e*=s;});
			auto res2 = thrust::inner_product(x.begin(), x.end(), y.begin(), T{}, thrust::plus<>{}, [] __device__(T const& t1, T const& t2) { return t1 * conj(t2); });
			BOOST_TEST(res == res2);
		}
		{
			BOOST_TEST( blas::dot(blas::C(x), x) == pow(blas::nrm2(x), 2.0) );
			BOOST_TEST( blas::dot(x, blas::C(x)) == pow(blas::nrm2(x), 2.0) );

			using blas::operators::operator, ;
			using blas::operators::operator*;
			using blas::operators::abs;
			using blas::operators::norm;
			using blas::operators::operator^;

			BOOST_TEST( (*x, x) == pow(abs(x), 2.0) );
			BOOST_TEST( (*x, x) == pow(abs(x), 2)   );
			BOOST_TEST( (*x, x) == norm(x)          );

			BOOST_TEST( (x, *x) == pow(abs(x), 2.0) );
			BOOST_TEST( (x, *x) == pow(abs(x), 2)   );
			BOOST_TEST( (x, *x) == norm(x)          );

			BOOST_TEST( (*x, x) == (x^2)            );
		}
	}
	{
		// T res;
		// blas::dot(blas::C(x), blas::C(y), res);
		multi::array<T, 2, Alloc> res({ 1, 1 }, 0.0);
		auto                      rr = blas::gemm(1.0, x.partitioned(1), blas::H(y.partitioned(1)), 0.0, res)[0][0];
		// {
		//  using blas::operators::operator,;
		//  auto res2 = +(x, blas::C(y));
		//  BOOST_TEST(res == res2);
		// }
		// {
		//  auto res2 = +blas::dot(x, blas::C(y));
		//  BOOST_TEST(res == res2);
		// }
		// {
		// //  auto [x2, A2, B2] = generate_ABx<complex, thrust::cuda::allocator<complex> >();
		// //  thrust::for_each(x2.begin(), x2.end(), [s] __device__ (T& e) {return e*=s;});
		//  auto res2 = thrust::inner_product(x.begin(), x.end(), y.begin(), T{}, thrust::plus<>{}, [] __device__ (T const& t1, T const& t2) {return t1*conj(t2);});
		//  BOOST_TEST(res == res2);
		// }
		// {
		//  BOOST_TEST( blas::dot(blas::C(x), x) == pow(blas::nrm2(x), 2.0) );
		//  BOOST_TEST( blas::dot(x, blas::C(x)) == pow(blas::nrm2(x), 2.0) );

		//  using blas::operators::operator,;
		//  using blas::operators::operator*;
		//  using blas::operators::abs;
		//  using blas::operators::norm;
		//  using blas::operators::operator^;

		//  BOOST_TEST( (*x, x) == pow(abs(x), 2.0) );
		//  BOOST_TEST( (*x, x) == pow(abs(x), 2)   );
		//  BOOST_TEST( (*x, x) == norm(x)          );

		//  BOOST_TEST( (x, *x) == pow(abs(x), 2.0) );
		//  BOOST_TEST( (x, *x) == pow(abs(x), 2)   );
		//  BOOST_TEST( (x, *x) == norm(x)          );

		//  BOOST_TEST( (*x, x) == (x^2)            );
		// }
	}
}

BOOST_AUTO_TEST_CASE(cublas_axpy_complex_one) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<complex, 1, Alloc> const x = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0},
		{4.1, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc> y = {
		{ 2.1, 0.0},
		{ 4.1, 0.0},
		{ 6.1, 0.0},
		{11.0, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming

	blas::axpy(1.0, x, y);
	std::cout << y[0] << std::endl;
	BOOST_TEST( static_cast<complex>(y[0]) == 3.2 + I*0.0 );
	{
		multi::array<complex, 1, Alloc> yy = {
			{ 2.1, 0.0},
			{ 4.1, 0.0},
			{ 6.1, 0.0},
			{11.0, 0.0}
		};

		thrust::transform(x.begin(), x.end(), yy.begin(), yy.begin(), [] __device__(auto const& ex, auto const& ey) { return ex + ey; });
		BOOST_TEST( yy == y );  // , boost::test_tools::per_element() );
	}
	{
		multi::array<complex, 1, Alloc> yy = {
			{ 2.1, 0.0},
			{ 4.1, 0.0},
			{ 6.1, 0.0},
			{11.0, 0.0}
		};
		using blas::operators::operator+=;
		yy += x;
		BOOST_TEST( yy == y );
	}
}

BOOST_AUTO_TEST_CASE(cublas_axpy_complex_mone) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<complex, 1, Alloc> const x = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0},
		{4.1, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc> y = {
		{ 2.1, 0.0},
		{ 4.1, 0.0},
		{ 6.1, 0.0},
		{11.0, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming

	blas::axpy(-1.0, x, y);
	std::cout << y[0] << std::endl;
	BOOST_TEST( static_cast<complex>(y[0]) == 1.0 + I*0.0 );
	{
		multi::array<complex, 1, Alloc> yy = {
			{ 2.1, 0.0},
			{ 4.1, 0.0},
			{ 6.1, 0.0},
			{11.0, 0.0}
		};  // NOLINT(readability-identifier-length) BLAS naming
		thrust::transform(x.begin(), x.end(), yy.begin(), yy.begin(), [] __host__ __device__(T ex, T ey) { return -1.0 * ex + ey; });
		BOOST_TEST( yy == y );  // boost::test_tools::per_element() );
	}
	{
		multi::array<complex, 1, Alloc> yy = {
			{ 2.1, 0.0},
			{ 4.1, 0.0},
			{ 6.1, 0.0},
			{11.0, 0.0}
		};
		using blas::operators::operator-=;
		yy -= x;
		BOOST_TEST( yy == y );
	}
	{
		multi::array<complex, 1, Alloc> yy = {
			{ 2.1, 0.0},
			{ 4.1, 0.0},
			{ 6.1, 0.0},
			{11.0, 0.0}
		};
		using blas::operators::operator-=;
		yy -= x;
		yy -= y;
		using blas::operators::norm;
		BOOST_TEST( norm(yy) == 0 );
		using blas::operators::operator==;
		BOOST_TEST( operator==(yy, 0) );
		BOOST_TEST( yy == 0 );
	}
}

BOOST_AUTO_TEST_CASE(cublas_axpy_complex_alpha) {
	namespace blas = multi::blas;
	complex const I{ 0.0, 1.0 };

	using T     = complex;
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<complex, 1, Alloc> const x = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0},
		{4.1, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc> y = {
		{ 2.1, 0.0},
		{ 4.1, 0.0},
		{ 6.1, 0.0},
		{11.0, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming

	blas::axpy(3.0, x, y);
	std::cout << y[0] << std::endl;
	BOOST_TEST( static_cast<complex>(y[0]) == 5.4 + I*0.0 );
	// {
	//  multi::array<complex, 1, Alloc> yy = {
	//      { 2.1, 0.0},
	//      { 4.1, 0.0},
	//      { 6.1, 0.0},
	//      {11.0, 0.0}
	//  };  // NOLINT(readability-identifier-length) BLAS naming
	//  thrust::transform(x.begin(), x.end(), yy.begin(), yy.begin(), [aa = 3.0] __device__(T ex, T ey) { return aa * ex + ey; });
	//  BOOST_TEST( yy == y , boost::test_tools::per_element() );
	// }
	{
		multi::array<complex, 1, Alloc> yy = {
			{ 2.1, 0.0},
			{ 4.1, 0.0},
			{ 6.1, 0.0},
			{11.0, 0.0}
		};
		using blas::operators::operator+=;
		using blas::operators::operator*;
		yy += 3.0 * x;
		BOOST_TEST( yy == y );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemv_conj_complex_zero) {
	namespace blas = multi::blas;
	using T        = complex;
	complex const I{ 0.0, 1.0 };
	using Alloc = thrust::cuda::allocator<complex>;

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{ { 9.0, 0.0 }, { 24.0, 0.0 }, { 30.0, 0.0 }, { 9.0, 0.0 }},
		{ { 4.0, 0.0 }, { 10.0, 0.0 }, { 12.0, 0.0 }, { 7.0, 0.0 }},
		{{ 14.0, 0.0 }, { 16.0, 0.0 }, { 36.0, 0.0 }, { 1.0, 0.0 }},
	};
	multi::array<complex, 1, Alloc> const x = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0},
		{4.1, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc> y = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0}
	};  // NOLINT(readability-identifier-length) BLAS naming
	blas::gemv(1.0, A, x, 0.0, y);
	{

		multi::array<complex, 1, Alloc> yy = {
			{1.1, 0.0},
			{2.1, 0.0},
			{3.1, 0.0}
		};  // NOLINT(readability-identifier-length) BLAS naming
		std::transform(begin(A), end(A), begin(yy), [&x](auto const& Ac) { return blas::dot(Ac, x); });

		BOOST_TEST( std::abs(static_cast<complex>(y[0]).real() - static_cast<complex>(yy[0]).real()) < 1e-10 );
		BOOST_TEST( std::abs(static_cast<complex>(y[1]).imag() - static_cast<complex>(yy[1]).imag()) < 1e-10 );
		BOOST_TEST( std::abs(static_cast<complex>(y[2]).real() - static_cast<complex>(yy[2]).real()) < 1e-10 );
	}
	{
		multi::array<complex, 1, Alloc> yy = {
			{1.1, 0.0},
			{2.1, 0.0},
			{3.1, 0.0}
		};  // NOLINT(readability-identifier-length) BLAS naming
		yy = blas::gemv(1.0, A, x);
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		multi::array<complex, 1, Alloc> yy = blas::gemv(1.0, A, x);
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		using blas::operators::operator%;

		multi::array<complex, 1, Alloc> yy = {
			{1.1, 0.0},
			{2.1, 0.0},
			{3.1, 0.0}
		};  // NOLINT(readability-identifier-length) BLAS naming
		yy = A % x;
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemv_complex_conj_zero) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{ 9.0 + I * 0.0, 24.0 + I * 0.0, 30.0 + I * 0.0, 9.0 + I * 0.0},
		{ 4.0 + I * 0.0, 10.0 + I * 0.0, 12.0 + I * 0.0, 7.0 + I * 0.0},
		{14.0 + I * 0.0, 16.0 + I * 0.0, 36.0 + I * 0.0, 1.0 + I * 0.0},
	};
	multi::array<complex, 1, Alloc> const x = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0 };                 // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc>       y = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
	blas::gemv(1.0, blas::T(A), x, 0.0, y);
	{
		multi::array<complex, 1, Alloc> yy = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
		using blas::operators::operator*;
	
		std::transform(begin(A.transposed()), end(A.transposed()), begin(yy), [&x](auto const& Ac) { return blas::dot(Ac, x); });

		BOOST_REQUIRE_CLOSE(static_cast<complex>(y[0]).real(), static_cast<complex>(yy[0]).real(), 1e-7);
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		multi::array<complex, 1, Alloc> yy = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
		yy                                 = blas::gemv(1.0, blas::T(A), x);
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		multi::array<complex, 1, Alloc> yy = blas::gemv(1.0, blas::T(A), x);
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		using blas::operators::operator%;

		multi::array<complex, 1, Alloc> yy = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
		yy                                 = ~A % x;
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemv_complex_zero) {
	namespace blas = multi::blas;
	using T        = complex;
	complex const I{ 0.0, 1.0 };
	using Alloc = thrust::cuda::allocator<complex>;

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{ { 9.0, 0.0 }, { 24.0, 0.0 }, { 30.0, 0.0 }, { 9.0, 0.0 }},
		{ { 4.0, 0.0 }, { 10.0, 0.0 }, { 12.0, 0.0 }, { 7.0, 0.0 }},
		{{ 14.0, 0.0 }, { 16.0, 0.0 }, { 36.0, 0.0 }, { 1.0, 0.0 }},
	};
	multi::array<complex, 1, Alloc> const x = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0},
		{4.1, 0.0}
	};

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc> y = {
		{1.1, 0.0},
		{2.1, 0.0},
		{3.1, 0.0}
	};
	blas::gemv(1.0, blas::J(A), x, 0.0, y);
	{
		multi::array<complex, 1, Alloc> yy = {
			{1.1, 0.0},
			{2.1, 0.0},
			{3.1, 0.0}
		};
		std::transform(begin(A), end(A), begin(yy), [&x](auto const& Ac) {
			using blas::operators::operator*;  // nvcc 11.8 needs this to be inside lambda
			return blas::dot(*Ac, x); }
		);

		BOOST_TEST( abs( static_cast<complex>(y[0]) - static_cast<complex>(yy[0])) < 1e-7 );
		BOOST_TEST( abs( static_cast<complex>(y[1]) - static_cast<complex>(yy[1])) < 1e-7 );
		BOOST_TEST( abs( static_cast<complex>(y[2]) - static_cast<complex>(yy[2])) < 1e-7 );
	}
	{
		multi::array<complex, 1, Alloc> yy = {
			{1.1, 0.0},
			{2.1, 0.0},
			{3.1, 0.0}
		};
		yy = blas::gemv(1.0, blas::J(A), x);
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		multi::array<complex, 1, Alloc> yy = blas::gemv(1.0, blas::J(A), x);
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		using blas::operators::operator%;
		using blas::operators::operator*;

		multi::array<complex, 1, Alloc> yy = {
			{1.1, 0.0},
			{2.1, 0.0},
			{3.1, 0.0}
		};  // NOLINT(readability-identifier-length) BLAS naming
		yy = *A % x;
		BOOST_TEST( static_cast<complex>(y[0]) == static_cast<complex>(yy[0]) );
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemv_complex_conjtrans_zero) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = std::allocator<complex>;  // thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{ 9.0 + I * 0.0, 24.0 + I * 0.0, 30.0 + I * 0.0, 9.0 + I * 0.0},
		{ 4.0 + I * 0.0, 10.0 + I * 0.0, 12.0 + I * 0.0, 7.0 + I * 0.0},
		{14.0 + I * 0.0, 16.0 + I * 0.0, 36.0 + I * 0.0, 1.0 + I * 0.0},
	};
	multi::array<complex, 1, Alloc> const x = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0 };                 // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc>       y = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming

	// blas::gemv(1.0, blas::H(A), x, 0.0, y);

	{
	// TODO(correaa) MKL gives an error here
	#if 0
		multi::array<complex, 1, Alloc> yy = { 1.1 + I* 0.0, 2.1 +I* 0.0, 3.1 + I* 0.0, 6.7 + I*0.0 };  // NOLINT(readability-identifier-length) BLAS naming
		std::transform(begin(transposed(A)), end(transposed(A)), begin(yy), [&x] (auto const& Ac) {
			using blas::operators::operator*;  // nvcc 11.8 needs this to be inside lambda
			return blas::dot(*Ac, x);}
		);

		BOOST_REQUIRE_CLOSE( static_cast<complex>(yy[0]).real() ,  61.7, 1.e-7  );
		BOOST_REQUIRE_CLOSE( static_cast<complex>(yy[1]).real() ,  97.0, 1.e-7  );
		BOOST_REQUIRE_CLOSE( static_cast<complex>(yy[2]).real() , 169.8, 1.e-7  );
		BOOST_REQUIRE_CLOSE( static_cast<complex>(yy[3]).real() ,  27.7, 1.e-7  );

		using blas::operators::operator*;
		BOOST_REQUIRE_CLOSE( static_cast<complex>(yy[0]).real() , (+blas::dot(*(~A)[0], x)).real() , 1.e-7  );
		BOOST_REQUIRE_CLOSE( static_cast<complex>(yy[1]).real() , (+blas::dot(*(~A)[1], x)).real() , 1.e-7  );
	#endif
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemv_complex_trans_one) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{ 9.0 + I * 0.0, 24.0 + I * 0.0, 30.0 + I * 0.0, 9.0 + I * 0.0},
		{ 4.0 + I * 0.0, 10.0 + I * 0.0, 12.0 + I * 0.0, 7.0 + I * 0.0},
		{14.0 + I * 0.0, 16.0 + I * 0.0, 36.0 + I * 0.0, 1.0 + I * 0.0},
	};
	multi::array<complex, 1, Alloc> const x = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0 };                 // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1, Alloc>       y = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
	blas::gemv(3.0 + I * 4.0, blas::T(A), x, 1.0, y);
	{
		multi::array<complex, 1, Alloc> yy = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
		// using blas::operators::operator*;
		std::transform(begin(transposed(A)), end(transposed(A)), begin(yy), begin(yy), [&x, aa = 3.0 + I * 4.0, bb = 1.0](auto const& Ac, complex e) { return aa * blas::dot(Ac, x) + bb * e; });

		BOOST_REQUIRE_CLOSE(static_cast<complex>(y[0]).real(), static_cast<complex>(yy[0]).real(), 1e-7);
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		multi::array<complex, 1, Alloc> yy = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
		yy += blas::gemv(3.0 + I * 4.0, blas::T(A), x);

		BOOST_REQUIRE_CLOSE(static_cast<complex>(y[0]).real(), static_cast<complex>(yy[0]).real(), 1e-7);
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
	{
		multi::array<complex, 1, Alloc> yy = { 1.1 + I * 0.0, 2.1 + I * 0.0, 3.1 + I * 0.0, 6.7 + I * 0.0 };  // NOLINT(readability-identifier-length) BLAS naming
		using blas::operators::operator*;
		yy += (3.0 + I * 4.0) * ~A % x;

		BOOST_REQUIRE_CLOSE(static_cast<complex>(y[0]).real(), static_cast<complex>(yy[0]).real(), 1e-7);
		BOOST_TEST( static_cast<complex>(y[1]) == static_cast<complex>(yy[1]) );
		BOOST_TEST( static_cast<complex>(y[2]) == static_cast<complex>(yy[2]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_trans_none) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		blas::gemm({ 1.0, 0.0 }, A, B, { 0.0, 0.0 }, C);

		// std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		//  [&A, aa=1.0, bb=0.0] (auto const& Bc, auto&& Cc) {return blas::gemv(aa, A, Bc, bb, std::move(Cc));}
		// );
		std::transform(begin(A), end(A), begin(C_copy), end(C_copy),
		               [&B, aa = 1.0, bb = 0](auto const& Ar, auto&& Cr) { return blas::gemv(aa, blas::T(B), Ar, bb, std::move(Cr)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C                                      = blas::gemm(1.0 + I * 0.0, A, B);

		// std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		//  [&A, aa=1.0, bb=0.0] (auto const& Bc, auto&& Cc) {return blas::gemv(aa, A, Bc, bb, std::move(Cc));}
		// );
		std::transform(begin(A), end(A), begin(C_copy), begin(C_copy), [&B, aa = 1.0, bb = 0.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, blas::T(B), Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C += blas::gemm(1.0 + I * 0.0, A, B);

		std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		               [&A, aa = 1.0, bb = 1.0](auto const& Bc, auto&& Cc) { return blas::gemv(aa, A, Bc, bb, std::move(Cc)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		using blas::operators::operator*;
		using blas::operators::operator+=;
		C += A * B;

		std::transform(begin(A), end(A), begin(C_copy), begin(C_copy), [&B, aa = 1.0, bb = 1.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, blas::T(B), Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_trans_second) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		blas::gemm({ 1.0, 0.0 }, A, blas::T(B), { 0.0, 0.0 }, C);

		std::transform(begin(B), end(B), begin(transposed(C_copy)), begin(transposed(C_copy)),
		               [&A, aa = 1.0, bb = 0.0](auto const& Bc, auto&& Cc) { return blas::gemv(aa, A, Bc, bb, std::move(Cc)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C                                      = blas::gemm(1.0 + I * 0.0, A, blas::T(B));

		// std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		//  [&A, aa=1.0, bb=0.0] (auto const& Bc, auto&& Cc) {return blas::gemv(aa, A, Bc, bb, std::move(Cc));}
		// );
		std::transform(begin(A), end(A), begin(C_copy), begin(C_copy), [&B, aa = 1.0, bb = 0.0](auto const& Ac, auto&& Cr) {
			return blas::gemv(aa, B, Ac, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C += blas::gemm(1.0 + I * 0.0, A, blas::T(B));

		std::transform(begin(B), end(B), begin(transposed(C_copy)), begin(transposed(C_copy)),
		               [&A, aa = 1.0, bb = 1.0](auto const& Bc, auto&& Cc) { return blas::gemv(aa, A, Bc, bb, std::move(Cc)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		using blas::operators::operator*;
		using blas::operators::operator+=;
		C += A * ~B;

		std::transform(begin(A), end(A), begin(C_copy), begin(C_copy), [&B, aa = 1.0, bb = 1.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, B, Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		using blas::operators::operator*;
		using blas::operators::operator+=;
		C += 2.0 * (A * ~B);

		std::transform(begin(A), end(A), begin(C_copy), begin(C_copy), [&B, aa = 2.0, bb = 1.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, B, Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_trans_first) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		blas::gemm({ 1.0, 0.0 }, blas::T(A), B, { 0.0, 0.0 }, C);

		std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		               [&A, aa = 1.0, bb = 0.0](auto const& Bc, auto&& Cc) { return blas::gemv(aa, blas::T(A), Bc, bb, std::move(Cc)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C                                      = blas::gemm(1.0 + I * 0.0, blas::T(A), B);

		// std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		//  [&A, aa=1.0, bb=0.0] (auto const& Bc, auto&& Cc) {return blas::gemv(aa, A, Bc, bb, std::move(Cc));}
		// );
		std::transform(begin(transposed(A)), end(transposed(A)), begin(C_copy), begin(C_copy), [&B, aa = 1.0, bb = 0.0](auto const& Ac, auto&& Cr) {
			return blas::gemv(aa, blas::T(B), Ac, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C += blas::gemm(1.0 + I * 0.0, blas::T(A), B);

		std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		               [&A, aa = 1.0, bb = 1.0](auto const& Bc, auto&& Cc) { return blas::gemv(aa, blas::T(A), Bc, bb, std::move(Cc)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		using blas::operators::operator*;
		using blas::operators::operator+=;
		C += ~A * B;

		std::transform(begin(transposed(A)), end(transposed(A)), begin(C_copy), begin(C_copy), [&B, aa = 1.0, bb = 1.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, blas::T(B), Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		using blas::operators::operator*;
		using blas::operators::operator+=;
		C += 2.0 * (~A * B);

		std::transform(begin(A.transposed()), end(A.transposed()), begin(C_copy), begin(C_copy), [&B, aa = 2.0, bb = 1.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, blas::T(B), Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_trans_both) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		blas::gemm({ 1.0, 0.0 }, blas::T(A), blas::T(B), { 0.0, 0.0 }, C);

		std::transform(B.begin(), B.end(), C_copy.transposed().begin(), C_copy.transposed().begin(),
		               [&A, aa = 1.0, bb = 0.0](auto const& Br, auto&& Cc) { return blas::gemv(aa, blas::T(A), Br, bb, std::move(Cc)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C                                      = blas::gemm(1.0 + I * 0.0, blas::T(A), blas::T(B));

		// std::transform(begin(transposed(B)), end(transposed(B)), begin(transposed(C_copy)), begin(transposed(C_copy)),
		//  [&A, aa=1.0, bb=0.0] (auto const& Bc, auto&& Cc) {return blas::gemv(aa, A, Bc, bb, std::move(Cc));}
		// );
		std::transform(A.transposed().begin(), A.transposed().end(), C_copy.begin(), C_copy.begin(), [&B, aa = 1.0, bb = 0.0](auto const& Ac, auto&& Cr) {
			return blas::gemv(aa, B, Ac, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		C += blas::gemm(1.0 + I * 0.0, blas::T(A), blas::T(B));

		std::transform(begin(B), end(B), C_copy.transposed().begin(), C_copy.transposed().begin(),
		               [&A, aa = 1.0, bb = 1.0](auto const& Br, auto&& Cc) { return blas::gemv(aa, blas::T(A), Br, bb, std::move(Cc)); }
		);

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		using blas::operators::operator*;
		using blas::operators::operator+=;
		C += ~A * ~B;

		std::transform(A.transposed().begin(), A.transposed().end(), C_copy.begin(), C_copy.begin(), [&B, aa = 1.0, bb = 1.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, B, Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            C_copy = C;
		using blas::operators::operator*;
		using blas::operators::operator+=;
		C += 2.0 * (~A * ~B);

		std::transform(A.transposed().begin(), A.transposed().end(), begin(C_copy), begin(C_copy), [&B, aa = 2.0, bb = 1.0](auto const& Ar, auto&& Cr) {
			return blas::gemv(aa, B, Ar, bb, std::move(Cr));
		});

		BOOST_TEST( static_cast<complex>(C_copy[1][0]) == static_cast<complex>(C[1][0]) );
		BOOST_TEST( static_cast<complex>(C_copy[0][1]) == static_cast<complex>(C[0][1]) );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_conj_second) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = std::allocator<complex>;  // thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            CC     = C;
		auto                            C_copy = CC;
		// blas::gemm({1.0, 0.0}, A, blas::J(B), {0.0, 0.0}, C);
		blas::gemm({ 1.0, 0.0 }, blas::T(B), blas::H(A), { 0.0, 0.0 }, C_copy);
		{
			auto const [is, js] = C.extensions();
			for(auto i : is) {
				for(auto j : js) {
					C[i][j] *= 0.0;
					for(auto k : B.extension()) {
						C[i][j] += A[i][k] * conj(B[k][j]);
					}
				}
			}
		}
	// TODO(correaa) MKL gives an error here
	// unknown location(0): fatal error: in "cublas_one_gemv_complex_conjtrans_zero": memory access violation at address: 0x00000007: no mapping at fault address
	#if 0
		{
			std::transform(begin(A), end(A), begin(CC), begin(CC), [BT = transposed(B)](auto const& Ar, auto&& Cr) {
				return std::transform(
					begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& Bc, auto&& Ce) {
						return 1.0*blas::dot(Ar, blas::C(Bc)) + 0.0*Ce;
					}
				), std::move(Cr);
			});
		}
		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );

		BOOST_TEST( static_cast<complex>(C_copy[1][0]).real() == +static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(C_copy[1][0]).imag() == -static_cast<complex>(C[0][1]).imag() );
	#endif
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_conj_first) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = std::allocator<complex>;  // thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            CC     = C;
		auto                            C_copy = CC;
		// blas::gemm({1.0, 0.0}, blas::J(A), B, {0.0, 0.0}, C);
		// blas::gemm({1.0, 0.0}, blas::T(B), blas::H(A), {0.0, 0.0}, C_copy);
		// {
		//  auto const [is, js] = C.extensions();
		//  for(auto i : is) {
		//      for(auto j : js) {
		//          C[i][j] *= 0.0;
		//          for(auto k : B.extension()) {
		//              C[i][j] += A[i][k]*conj(B[k][j]);
		//          }
		//      }
		//  }
		// }
		// {
		//  std::transform(begin(A), end(A), begin(CC), begin(CC), [BT = transposed(B)](auto const& Ar, auto&& Cr) {
		//      return std::transform(
		//          begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& BCr, auto&& Ce) {
		//              return 1.0*blas::dot(Ar, blas::C(BCr)) + 0.0*Ce;
		//          }
		//      ), std::move(Cr);
		//  });
		// }
		// BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		// BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		// BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		// BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );

		// BOOST_TEST( static_cast<complex>(C_copy[1][0]).real() == +static_cast<complex>(C[0][1]).real() );
		// BOOST_TEST( static_cast<complex>(C_copy[1][0]).imag() == -static_cast<complex>(C[0][1]).imag() );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_conj_both) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = std::allocator<complex>;  // thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	{
		multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
		auto                            CC     = C;
		auto                            C_copy = CC;
		//  blas::gemm({1.0, 0.0}, blas::J(A), blas::J(B), {0.0, 0.0}, C);
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_herm_second) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
	blas::gemm({ 1.0, 0.0 }, A, blas::H(B), { 0.0, 0.0 }, C);
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });

		std::transform(
			begin(A), end(A), begin(CC), begin(CC),
			[&B, aa = 1.0, bb = 0.0](auto const& Ar, auto&& Cr) {
				return blas::gemv(aa, blas::J(B), Ar, bb, std::move(Cr));
			}
		);

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		CC = blas::gemm({ 1.0, 0.0 }, A, blas::H(B));

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		using blas::operators::operator*;
		using blas::operators::operator~;
		CC = A * ~*B;

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_herm_second_plus) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
	blas::gemm({ 1.0, 0.0 }, A, blas::H(B), { 1.0, 0.0 }, C);
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });

		std::transform(
			begin(A), end(A), begin(CC), begin(CC),
			[&B, aa = 1.0, bb = 1.0](auto const& Ar, auto&& Cr) {
				return blas::gemv(aa, blas::J(B), Ar, bb, std::move(Cr));
			}
		);

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		CC += blas::gemm({ 1.0, 0.0 }, A, blas::H(B));

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		using blas::operators::operator*;
		using blas::operators::operator~;
		CC += A * ~*B;

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_herm_first) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
	blas::gemm({ 1.0, 0.0 }, blas::H(A), B, { 0.0, 0.0 }, C);
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//  auto const [is, js] = CC.extensions();
	//  for(auto i : is) {
	//      for(auto j : js) {
	//          CC[i][j] = 0.0;
	//          for(auto k : A.extension()) {
	//              CC[i][j] += 1.0*conj(A[k][i])*B[k][j] ;
	//          }
	//      }
	//  }
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });

		std::transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [BT = transposed(B)](auto const& Ac, auto&& Cr) {
			std::transform(begin(BT), end(BT), begin(Cr), begin(Cr), [&Ac](auto const& Bc, auto&& c) {
				return blas::dot(blas::C(Ac), Bc, std::move(c));
			});
			return std::move(Cr);
		});
		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		CC = blas::gemm({ 1.0, 0.0 }, blas::H(A), B);

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		using blas::operators::operator*;
		using blas::operators::operator~;
		CC = ~*A * B;

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_herm_both) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
	blas::gemm({ 1.0, 0.0 }, blas::H(A), blas::H(B), { 0.0, 0.0 }, C);
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//  auto const [is, js] = CC.extensions();
	//  for(auto i : is) {
	//      for(auto j : js) {
	//          CC[i][j] = 0.0;
	//          for(auto k : A.extension()) {
	//              CC[i][j] += 1.0*conj(A[k][i])*conj(B[j][k]) ;
	//          }
	//      }
	//  }
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  thrust::transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [BP = &B] __device__ (multi::array<complex, 2, Alloc>::reference const& Ac, multi::array<complex, 2, Alloc>::reference&& Cr) {
	//      thrust::transform(begin(*BP), end(*BP), begin(Cr), begin(Cr), [APc = &Ac] __device__  (multi::array<complex, 2, Alloc>::reference const& Bc, complex&& c) {
	//          return conj(thrust::inner_product(begin(*APc), end(*APc), begin(Bc), 0.0*c, std::plus<>{}, [] __device__ (complex const& a, complex const& b) {return a*b;}));
	//      //  return conj(+blas::dot(Ac, Bc, std::move(c)));
	//      });
	//      return std::move(Cr);
	//  });
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });

		CC = blas::gemm({ 1.0, 0.0 }, blas::H(A), blas::H(B));

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		using blas::operators::operator*;
		using blas::operators::operator~;
		CC = ~*A * ~*B;

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_trans_herm) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
	blas::gemm({ 1.0, 0.0 }, blas::T(A), blas::H(B), { 0.0, 0.0 }, C);
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//  auto const [is, js] = CC.extensions();
	//  for(auto i : is) {
	//      for(auto j : js) {
	//          CC[i][j] = 0.0;
	//          for(auto k : A.extension()) {
	//              CC[i][j] += 1.0*conj(A[k][i])*conj(B[j][k]) ;
	//          }
	//      }
	//  }
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  thrust::transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [BP = &B] __device__ (multi::array<complex, 2, Alloc>::reference const& Ac, multi::array<complex, 2, Alloc>::reference&& Cr) {
	//      thrust::transform(begin(*BP), end(*BP), begin(Cr), begin(Cr), [APc = &Ac] __device__  (multi::array<complex, 2, Alloc>::reference const& Bc, complex&& c) {
	//          return conj(thrust::inner_product(begin(*APc), end(*APc), begin(Bc), 0.0*c, std::plus<>{}, [] __device__ (complex const& a, complex const& b) {return a*b;}));
	//      //  return conj(+blas::dot(Ac, Bc, std::move(c)));
	//      });
	//      return std::move(Cr);
	//  });
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });

		CC = blas::gemm({ 1.0, 0.0 }, blas::T(A), blas::H(B));

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
	{
		multi::array<complex, 2, Alloc> CC({ 2, 2 }, { 3.0, 0.0 });
		using blas::operators::operator*;
		using blas::operators::operator~;
		CC = ~A * ~*B;

		BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
		BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

		BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
		BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	}
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_herm_trans) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming

	// blas::gemm({1.0, 0.0}, blas::H(A), blas::T(B), {0.0, 0.0}, C);
	//  {
	//   multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//   auto const [is, js] = CC.extensions();
	//   for(auto i : is) {
	//       for(auto j : js) {
	//           CC[i][j] = 0.0;
	//           for(auto k : A.extension()) {
	//               CC[i][j] += 1.0*conj(A[k][i])*conj(B[j][k]) ;
	//           }
	//       }
	//   }
	//   BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//   BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  thrust::transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [BP = &B] __device__ (multi::array<complex, 2, Alloc>::reference const& Ac, multi::array<complex, 2, Alloc>::reference&& Cr) {
	//      thrust::transform(begin(*BP), end(*BP), begin(Cr), begin(Cr), [APc = &Ac] __device__  (multi::array<complex, 2, Alloc>::reference const& Bc, complex&& c) {
	//          return conj(thrust::inner_product(begin(*APc), end(*APc), begin(Bc), 0.0*c, std::plus<>{}, [] __device__ (complex const& a, complex const& b) {return a*b;}));
	//      //  return conj(+blas::dot(Ac, Bc, std::move(c)));
	//      });
	//      return std::move(Cr);
	//  });
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  CC = blas::gemm({1.0, 0.0}, blas::H(A), blas::T(B));

	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//  using blas::operators::operator*;
	//  using blas::operators::operator~;
	//  CC = ~*A * ~B;

	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_conj_herm) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming

	// blas::gemm({1.0, 0.0}, blas::J(A), blas::H(B), {0.0, 0.0}, C);
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//  auto const [is, js] = CC.extensions();
	//  for(auto i : is) {
	//      for(auto j : js) {
	//          CC[i][j] = 0.0;
	//          for(auto k : A.extension()) {
	//              CC[i][j] += 1.0*conj(A[k][i])*conj(B[j][k]) ;
	//          }
	//      }
	//  }
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  thrust::transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [BP = &B] __device__ (multi::array<complex, 2, Alloc>::reference const& Ac, multi::array<complex, 2, Alloc>::reference&& Cr) {
	//      thrust::transform(begin(*BP), end(*BP), begin(Cr), begin(Cr), [APc = &Ac] __device__  (multi::array<complex, 2, Alloc>::reference const& Bc, complex&& c) {
	//          return conj(thrust::inner_product(begin(*APc), end(*APc), begin(Bc), 0.0*c, std::plus<>{}, [] __device__ (complex const& a, complex const& b) {return a*b;}));
	//      //  return conj(+blas::dot(Ac, Bc, std::move(c)));
	//      });
	//      return std::move(Cr);
	//  });
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  CC = blas::gemm({1.0, 0.0}, blas::T(A), blas::H(B));

	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//  using blas::operators::operator*;
	//  using blas::operators::operator~;
	//  CC = ~A * ~*B;

	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
}

BOOST_AUTO_TEST_CASE(cublas_one_gemm_complex_herm_conj) {
	namespace blas = multi::blas;
	using T        = complex;
	using Alloc    = thrust::cuda::allocator<complex>;
	complex const I{ 0.0, 1.0 };

	// NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 - 2.0 * I, 9.0 - 1.0 * I},
		{2.0 + 3.0 * I, 1.0 - 2.0 * I},
	};

	multi::array<complex, 2, Alloc> const B = {
		{3.0 - 4.0 * I, 19.0 - 1.0 * I},
		{1.0 + 5.0 * I,  8.0 - 8.0 * I},
	};
	multi::array<complex, 2, Alloc> C({ 2, 2 }, { 3.0, 0.0 });  // NOLINT(readability-identifier-length) conventional BLAS naming
	                                                            // blas::gemm({1.0, 0.0}, blas::H(A), blas::J(B), {0.0, 0.0}, C);
	                                                            //  {
	                                                            //   multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	                                                            //   auto const [is, js] = CC.extensions();
	                                                            //   for(auto i : is) {
	                                                            //       for(auto j : js) {
	                                                            //           CC[i][j] = 0.0;
	                                                            //           for(auto k : A.extension()) {
	                                                            //               CC[i][j] += 1.0*conj(A[k][i])*conj(B[j][k]) ;
	                                                            //           }
	                                                            //       }
	                                                            //   }
	                                                            //   BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	                                                            //   BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  thrust::transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [BP = &B] __device__ (multi::array<complex, 2, Alloc>::reference const& Ac, multi::array<complex, 2, Alloc>::reference&& Cr) {
	//      thrust::transform(begin(*BP), end(*BP), begin(Cr), begin(Cr), [APc = &Ac] __device__  (multi::array<complex, 2, Alloc>::reference const& Bc, complex&& c) {
	//          return conj(thrust::inner_product(begin(*APc), end(*APc), begin(Bc), 0.0*c, std::plus<>{}, [] __device__ (complex const& a, complex const& b) {return a*b;}));
	//      //  return conj(+blas::dot(Ac, Bc, std::move(c)));
	//      });
	//      return std::move(Cr);
	//  });
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});

	//  CC = blas::gemm({1.0, 0.0}, blas::T(A), blas::H(B));

	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
	// {
	//  multi::array<complex, 2, Alloc> CC({2, 2}, {3.0, 0.0});
	//  using blas::operators::operator*;
	//  using blas::operators::operator~;
	//  CC = ~A * ~*B;

	//  BOOST_TEST( static_cast<complex>(CC[1][0]).real() == static_cast<complex>(C[1][0]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[1][0]).imag() == static_cast<complex>(C[1][0]).imag() );

	//  BOOST_TEST( static_cast<complex>(CC[0][1]).real() == static_cast<complex>(C[0][1]).real() );
	//  BOOST_TEST( static_cast<complex>(CC[0][1]).imag() == static_cast<complex>(C[0][1]).imag() );
	// }
}

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<double>;
	complex const I{ 0.0, 1.0 };  // NOLINT(readability-identifier-length) imag unit
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<complex, 2, Alloc> const A = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 4.0 * I, 3.0 + 0.0 * I, 4.0 - 10.0 * I},
		{0.0 + 0.0 * I, 7.0 - 3.0 * I,  1.0 + 0.0 * I},
		{0.0 + 0.0 * I, 0.0 + 0.0 * I,  8.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> B = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 1.0 * I, 2.0 + 1.0 * I, 3.0 + 1.0 * I},
		{5.0 + 3.0 * I, 9.0 + 3.0 * I, 1.0 - 1.0 * I},
	};

	using multi::blas::trsm;

	blas::trsm(blas::side::left, { 1.0, 0.0 }, blas::U(A), blas::H(B));  // B†←A⁻¹.B†, B←B.A⁻¹†, B←(A⁻¹.B†)†
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[1][2]).imag(), -0.147059, 0.001);
}

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_no_const_UTH) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<double>;
	complex const I{ 0.0, 1.0 };  // NOLINT(readability-identifier-length) imag unit
	using Alloc = thrust::cuda::allocator<complex>;

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{ 1.0 + 4.0 * I, 0.0 + 0.0 * I, 0.0 - 0.0 * I},
		{ 3.0 + 0.0 * I, 7.0 - 3.0 * I, 0.0 + 0.0 * I},
		{4.0 - 10.0 * I, 1.0 + 0.0 * I, 8.0 - 2.0 * I},
	};

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> B = {
		{1.0 + 1.0 * I, 2.0 + 1.0 * I},
		{5.0 + 3.0 * I, 9.0 + 3.0 * I},
		{3.0 + 1.0 * I, 1.0 - 1.0 * I},
	};

	using multi::blas::trsm;

	blas::trsm(blas::side::left, { 1.0, 0.0 }, blas::U(blas::H(A)), B);
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[1][1]).imag(), -0.0811359, 0.001);
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[2][1]).imag(), -0.147059, 0.001);
}

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_gemm_check_no_const) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<double>;
	complex const I{ 0.0, 1.0 };  // NOLINT(readability-identifier-length) imag unit
	using Alloc = thrust::cuda::allocator<complex>;

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 + 4.0 * I, 3.0 + 0.0 * I, 4.0 - 10.0 * I},
		{0.0 + 0.0 * I, 7.0 - 3.0 * I,  1.0 + 0.0 * I},
		{0.0 + 0.0 * I, 0.0 + 0.0 * I,  8.0 - 2.0 * I},
	};

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> B = {
		{1.0 + 1.0 * I, 2.0 + 1.0 * I},
		{5.0 + 3.0 * I, 9.0 + 3.0 * I},
		{3.0 + 1.0 * I, 1.0 - 1.0 * I},
	};

	using multi::blas::trsm;

	blas::trsm(blas::side::left, { 1.0, 0.0 }, blas::U(A), B);  // B←A⁻¹.B, B†←A⁻¹.B†
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[2][1]).imag(), -0.0882353, 0.001);
}

// BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_gemm_check_no_const_conj_second) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
//  using Alloc = thrust::cuda::allocator<complex>;

//  // NOLINTNEXTLINE(readability-identifier-length) BLAS naming
//  multi::array<complex, 2, Alloc> const A = {
//      { 1.0 + 4.0*I, 3.0 + 0.0*I,  4.0 - 10.0*I},
//      { 0.0 + 0.0*I, 7.0 - 3.0*I,  1.0 +  0.0*I},
//      { 0.0 + 0.0*I, 0.0 + 0.0*I,  8.0 -  2.0*I},
//  };

//  // NOLINTNEXTLINE(readability-identifier-length) BLAS naming
//  multi::array<complex, 2, Alloc> B = {
//      {1.0 + 1.0*I, 2.0 + 1.0*I},
//      {5.0 + 3.0*I, 9.0 + 3.0*I},
//      {3.0 + 1.0*I, 1.0 - 1.0*I},
//  };

//  using multi::blas::trsm;

//  blas::trsm(blas::side::left, {1.0, 0.0}, blas::U(A), blas::J(B));  // B*←A⁻¹.B*, B^T←A⁻¹.B^T
//  BOOST_REQUIRE_CLOSE( static_cast<complex>(B[2][1]).imag() , -0.0882353, 0.001);
// }

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_gemm_check_no_const_operator) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<double>;
	complex const I{ 0.0, 1.0 };  // NOLINT(readability-identifier-length) imag unit
	using Alloc = thrust::cuda::universal_allocator<complex>;

	multi::array<complex, 2> const A = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 4.0 * I, 3.0 + 0.0 * I, 4.0 - 10.0 * I},
		{0.0 + 0.0 * I, 7.0 - 3.0 * I,  1.0 + 0.0 * I},
		{0.0 + 0.0 * I, 0.0 + 0.0 * I,  8.0 - 2.0 * I},
	};
	multi::array<complex, 2> B = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 1.0 * I, 2.0 + 1.0 * I},
		{5.0 + 3.0 * I, 9.0 + 3.0 * I},
		{3.0 + 1.0 * I, 1.0 - 1.0 * I},
	};

	using blas::operators::operator|=;
	using blas::operators::U;
	B |= U(A);  // B←A⁻¹.B, B†←A⁻¹.B†
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[2][1]).imag(), -0.0882353, 0.001);
}

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_gemm_check_no_const_right) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<double>;
	complex const I{ 0.0, 1.0 };  // NOLINT(readability-identifier-length) imag unit
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<complex, 2, Alloc> const A = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 4.0 * I, 3.0 + 0.0 * I, 4.0 - 10.0 * I},
		{0.0 + 0.0 * I, 7.0 - 3.0 * I,  1.0 + 0.0 * I},
		{0.0 + 0.0 * I, 0.0 + 0.0 * I,  8.0 - 2.0 * I},
	};

	multi::array<complex, 2, Alloc> B = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 1.0 * I, 2.0 + 1.0 * I, 3.0 + 1.0 * I},
		{5.0 + 3.0 * I, 9.0 + 3.0 * I, 1.0 - 1.0 * I},
	};

	using multi::blas::trsm;

	blas::trsm(blas::side::right, { 1.0, 0.0 }, blas::U(A), B);  // B←B.A⁻¹, B←B/A, B†←A⁻¹†.B†
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[1][2]).imag(), 1.60142, 0.001);
}

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_gemm_check_no_const_right_LT) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<double>;
	complex const I{ 0.0, 1.0 };  // NOLINT(readability-identifier-length) imag unit
	using Alloc = thrust::cuda::allocator<complex>;

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> const A = {
		{1.0 + 4.0 * I, 3.0 + 0.0 * I, 4.0 - 10.0 * I},
		{0.0 + 0.0 * I, 7.0 - 3.0 * I,  1.0 + 0.0 * I},
		{0.0 + 0.0 * I, 0.0 + 0.0 * I,  8.0 - 2.0 * I},
	};

	// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
	multi::array<complex, 2, Alloc> B = {
		{1.0 + 1.0 * I, 2.0 + 1.0 * I, 3.0 + 1.0 * I},
		{5.0 + 3.0 * I, 9.0 + 3.0 * I, 1.0 - 1.0 * I},
	};

	using multi::blas::trsm;

	blas::trsm(blas::side::right, { 1.0, 0.0 }, blas::L(blas::T(A)), B);  // B←B.Aᵀ⁻¹, B←B/Aᵀ, B†←Aᵀ⁻¹†.B†, Bᵀ←A⁻¹.Bᵀ, Bᵀ←Bᵀ\A
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[1][2]).imag(), -0.0882353, 0.001);
}

// BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_gemm_check_no_const_right_LH) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>; complex const I{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit
//  using Alloc = thrust::cuda::allocator<complex>;

//  // NOLINTNEXTLINE(readability-identifier-length) BLAS naming
//  multi::array<complex, 2, Alloc> const A = {
//      { 1.0 + 4.0*I, 3.0 + 0.0*I,  4.0 - 10.0*I},
//      { 0.0 + 0.0*I, 7.0 - 3.0*I,  1.0 +  0.0*I},
//      { 0.0 + 0.0*I, 0.0 + 0.0*I,  8.0 -  2.0*I},
//  };

//  // NOLINTNEXTLINE(readability-identifier-length) BLAS naming
//  multi::array<complex, 2, Alloc> B = {
//      { 1.0 + 1.0*I, 2.0 + 1.0*I, 3.0 + 1.0*I},
//      { 5.0 + 3.0*I, 9.0 + 3.0*I, 1.0 - 1.0*I},
//  };

//  using multi::blas::trsm;

//  blas::trsm(blas::side::right, {1.0, 0.0}, blas::U(blas::J(A)), B);  // B←B.A*⁻¹, B←B/A*, B*←B*.A⁻¹
//  BOOST_REQUIRE_CLOSE( static_cast<complex>(B[1][2]).imag(), -0.0882353, 0.001);
// }

BOOST_AUTO_TEST_CASE(UTA_blas_trsm_complex_nonsquare_default_diagonal_gemm_check_no_const_right_operator) {
	namespace blas = multi::blas;
	using complex  = thrust::complex<double>;
	complex const I{ 0.0, 1.0 };  // NOLINT(readability-identifier-length) imag unit
	using Alloc = thrust::cuda::allocator<complex>;

	multi::array<complex, 2, Alloc> const A = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 4.0 * I, 3.0 + 0.0 * I, 4.0 - 10.0 * I},
		{0.0 + 0.0 * I, 7.0 - 3.0 * I,  1.0 + 0.0 * I},
		{0.0 + 0.0 * I, 0.0 + 0.0 * I,  8.0 - 2.0 * I},
	};
	multi::array<complex, 2, Alloc> B = {
  // NOLINT(readability-identifier-length) BLAS naming
		{1.0 + 1.0 * I, 2.0 + 1.0 * I, 3.0 + 1.0 * I},
		{5.0 + 3.0 * I, 9.0 + 3.0 * I, 1.0 - 1.0 * I},
	};

	using multi::blas::trsm;

	using blas::operators::operator/=;
	B /= blas::U(A);
	BOOST_REQUIRE_CLOSE(static_cast<complex>(B[1][2]).imag(), 1.60142, 0.001);
}

#endif

return boost::report_errors();}

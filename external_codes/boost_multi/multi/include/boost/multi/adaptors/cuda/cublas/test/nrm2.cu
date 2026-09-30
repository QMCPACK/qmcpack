// Copyright 2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/cuda/cublas.hpp>

#include <boost/multi/adaptors/blas/dot.hpp>
#include <boost/multi/adaptors/blas/axpy.hpp>
#include <boost/multi/adaptors/blas/gemm.hpp>
#include <boost/multi/adaptors/blas/nrm2.hpp>
#include <boost/multi/adaptors/blas/scal.hpp>

#include <boost/multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

#include<numeric>

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

int main() {
// BOOST_AUTO_TEST_CASE(cublas_dot_out_param_complex_C) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>;
//  complex const I{0.0, 1.0};

//  multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
//  multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

//  complex res{0.0, 0.0};
//  blas::dot(blas::C(x), y, res);
// //  BOOST_TEST( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
// }

BOOST_AUTO_TEST_CASE(cublas_dot_out_array0D_complex_C) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>;
	complex const I{0.0, 1.0};

	multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

	multi::thrust::cuda::array<complex, 0> res{complex{0.0, 0.0}};
	blas::dot(blas::C(x), y, res);

	{
		multi::array<complex, 0> res_copy{complex{0.0, 0.0}};
		res_copy = res;
		BOOST_TEST(( *res_copy.base() == complex{14.0, 4.0} ));
	}
	{
		multi::array<complex, 0> res_copy{res};
		BOOST_TEST(( *res_copy.base() == complex{14.0, 4.0} ));
	}
}

BOOST_AUTO_TEST_CASE(cublas_dot_out_array0D_complex_C) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>;
	complex const I{0.0, 1.0};

	multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

	// multi::thrust::cuda::array<complex, 0> res{complex{0.0, 0.0}};

	{
		double res{0.0};
		multi::array<double, 1> const xx = {1.0, 2.0, 3.0};

		res = blas::nrm2(xx);
		BOOST_TEST( std::abs( res - std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) ) < 1e-7 );
	}
	{
		double res{0.0};
		multi::array<double, 1> const xx = {1.0, 2.0, 3.0};

		blas::nrm2(xx, res);
		BOOST_TEST( std::abs( res - std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) ) < 1e-7 );
	}
	{
		multi::array<double, 0> res{0.0};
		multi::array<double, 1> const xx = {1.0, 2.0, 3.0};

		blas::nrm2(xx, res);
		BOOST_TEST( std::abs( *res.base() - std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) ) < 1e-7 );
	}
	{
		multi::array<double, 0> res{0.0};
		multi::array<std::complex<double>, 1> const xx = {std::complex<double>{1.0, 2.0}, std::complex<double>{3.0, 4.0}, std::complex<double>{5.0, 6.0}};

		blas::nrm2(xx, res);
		BOOST_TEST( std::abs( *res.base() - std::sqrt(std::norm(std::complex<double>{1.0, 2.0}) + std::norm(std::complex<double>{3.0, 4.0}) + std::norm(std::complex<double>{5.0, 6.0})) ) < 1e-7 );
	}
	{
		multi::thrust::cuda::array<double, 0> res{0.0};
		multi::thrust::cuda::array<double, 1> const xx = {1.0, 2.0, 3.0};

		blas::nrm2(xx, res);
		BOOST_TEST( *res.base() == std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) );
	}
	{
		multi::thrust::cuda::array<double, 0> res{0.0};
		multi::thrust::cuda::array<double, 1> const xx = {1.0, 2.0, 3.0};

		res = blas::nrm2(xx);
		BOOST_TEST( *res.base() == std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) );

		double res2 = blas::nrm2(xx);
		BOOST_TEST( res2 == std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) );
	}
	{
		multi::thrust::cuda::array<double, 0> res{0.0};
		multi::thrust::cuda::array<thrust::complex<double>, 1> const xx = {thrust::complex<double>(1.0, 2.0), thrust::complex<double>(3.0, 4.0), thrust::complex<double>(4.0, 5.0)};

		res = blas::nrm2(xx);
		BOOST_TEST( *res.base() == std::sqrt(std::norm(std::complex<double>(1.0, 2.0)) + std::norm(std::complex<double>(3.0, 4.0)) + std::norm(std::complex<double>(4.0, 5.0))) );

		double res2 = blas::nrm2(xx);
		BOOST_TEST( res2 == *res.base() );
	}
	{
		multi::thrust::cuda::array<double, 0> res{0.0};
		multi::thrust::cuda::array<thrust::complex<double>, 1> const xx = {thrust::complex<double>(1.0, 2.0), thrust::complex<double>(3.0, 4.0), thrust::complex<double>(4.0, 5.0)};

		blas::nrm2_n(xx.begin(), xx.size(), res.base());
		BOOST_TEST( *res.base() == std::sqrt(std::norm(std::complex<double>(1.0, 2.0)) + std::norm(std::complex<double>(3.0, 4.0)) + std::norm(std::complex<double>(4.0, 5.0))) );
	}
	{
		multi::thrust::cuda::array<double, 0> res{0.0};
		multi::thrust::cuda::array<double, 1> const xx = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

		blas::nrm2_n(xx.begin(), 3, res.base());
		BOOST_TEST( *res.base() == std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) );
	}

	// {
	//  multi::thrust::cuda::array<double, 0> res{0.0};
	//  multi::thrust::cuda::array<double, 1> const xx = {1.0, 2.0, 3.0};

	//  blas::nrm2(xx, res);
	//  BOOST_TEST( *res.base() == std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) );
	// }

	// multi::thrust::cuda::array<complex, 0> res{complex{0.0, 0.0}};
	// blas::dot(blas::C(x), y, res);

	// {
	//  multi::array<complex, 0> res_copy{complex{0.0, 0.0}};
	//  res_copy = res;
	//  BOOST_TEST(( *res_copy.base() == complex{14.0, 4.0} ));
	// }
	// {
	//  multi::array<complex, 0> res_copy{res};
	//  BOOST_TEST(( *res_copy.base() == complex{14.0, 4.0} ));
	// }
}

// BOOST_AUTO_TEST_CASE(blas_dot_functional_complex_C) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>;
//  complex const I{0.0, 1.0};

//  multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
//  multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

//  complex res = blas::dot(blas::C(x), y);
//  BOOST_TEST( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
// }

// BOOST_AUTO_TEST_CASE(blas_dot_functional_mutate_complex_C) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>;
//  complex const I{0.0, 1.0};

//  multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
//  multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

//  complex res;
//  res = blas::dot(blas::C(x), y);
//  BOOST_TEST( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
// }

return boost::report_errors();}

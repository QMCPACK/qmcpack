// Copyright 2020-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #include <boost/test/unit_test.hpp>  // TODO(correaa) convert into lightweight test

#include <boost/multi/adaptors/blas.hpp>
#include <boost/multi/adaptors/cuda.hpp>
#include <boost/multi/adaptors/cuda/cublas.hpp>
#include <boost/multi/array.hpp>

#include <chrono>
#include <cmath>  // for abs  // IWYU pragma: keep
// IWYU pragma: no_include <cstdlib>                                   // for abs
#include <exception>  // for exception
#include <iostream>
#include <random>
#include <string>
#include <utility>  // move

namespace multi = boost::multi;

dasdasdas

namespace {
class auto_timer {
	std::string label_;

	std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();

 public:
	explicit auto_timer(std::string label = {}) : label_{std::move(label)} {}

	auto_timer(auto_timer const&)                    = delete;
	auto_timer(auto_timer&&)                         = delete;

	auto operator=(auto_timer const&) -> auto_timer& = delete;
	auto operator=(auto_timer&&) -> auto_timer&      = delete;

	~auto_timer() {
		std::cerr << label_ << std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count() << " s (wall)\n";
	}
};
}  // namespace

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_gemm_complex_3x2_3x2 const) {
	using complex = std::complex<double>;
	complex const I{0, 1};
	namespace blas                   = multi::blas;
	multi::array<complex, 2> const a = {
		{1.0 + 2.0 * I, 5.0 + 2.0 * I},
		{9.0 - 1.0 * I, 9.0 + 1.0 * I},
		{1.0 + 1.0 * I, 2.0 + 2.0 * I}
	};
	multi::array<complex, 2> const b = {
		{11.0 - 2.0 * I, 5.0 + 2.0 * I},
		{ 7.0 - 3.0 * I, 2.0 + 1.0 * I},
		{ 8.0 - 1.0 * I, 1.0 + 1.0 * I}
	};
	{
		{
			multi::array<complex, 2> c({2, 2});
			c = blas::gemm(1., blas::H(a), b);  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE(c[1][0] == 125.0 - 84.0 * I);
		}
	}
	{
		multi::cuda::array<complex, 2> const a_gpu = a;
		multi::cuda::array<complex, 2> const b_gpu = b;
		{
			multi::cuda::array<complex, 2> c_gpu({2, 2});
			c_gpu = blas::gemm(1., blas::H(a_gpu), b_gpu);  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE(c_gpu[1][0] == 125.0 - 84.0 * I);
		}
		{
			auto c_gpu = +blas::gemm(1.0, blas::H(a_gpu), b_gpu);
			BOOST_REQUIRE(c_gpu[1][0] == 125.0 - 84.0 * I);
		}
	}
	{
		multi::cuda::managed::array<complex, 2> const a_gpu = a;
		multi::cuda::managed::array<complex, 2> const b_gpu = b;
		{
			multi::cuda::managed::array<complex, 2> c_gpu({2, 2});
			blas::gemm(1., blas::H(a_gpu), b_gpu, 0., c_gpu);  // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE(c_gpu[1][0] == 125.0 - 84.0 * I);
		}
		{
			auto c_gpu = +blas::gemm(1.0, blas::H(a_gpu), b_gpu);
			BOOST_REQUIRE(c_gpu[1][0] == 125.0 - 84.0 * I);
		}
	}
}

// BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_gemm_complex_3x2_3x2_with_context){
//   using complex = std::complex<double>; complex const I{0, 1};
//   namespace blas = multi::blas;
//   multi::array<complex, 2> const a = {
//       {1. + 2.*I, 5. + 2.*I},
//       {9. - 1.*I, 9. + 1.*I},
//       {1. + 1.*I, 2. + 2.*I}
//   };
//   multi::array<complex, 2> const b = {
//       { 11. - 2.*I, 5. + 2.*I},
//       {  7. - 3.*I, 2. + 1.*I},
//       {  8. - 1.*I, 1. + 1.*I}
//   };
//   {
//       {
//           multi::blas::context ctx;
//           multi::array<complex, 2> c({2, 2});
//           blas::gemm(ctx, 1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
//           BOOST_REQUIRE( c[1][0] == 125.-84.*I );
//       }
//   }
//   {
//       multi::cublas::context ctx;
//       multi::cuda::array<complex, 2> const a_gpu = a;
//       multi::cuda::array<complex, 2> const b_gpu = b;
//       {
//           multi::cuda::array<complex, 2> c_gpu({2, 2});
//           blas::gemm(ctx, 1., blas::H(a_gpu), b_gpu, 0., c_gpu); // c=ab, c⸆=b⸆a⸆
//           BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//       }
//       {
//           auto c_gpu =+ blas::gemm(&ctx, blas::H(a_gpu), b_gpu);
//           BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//       }
//   }
//   {
//       multi::cublas::context ctx;
//       multi::cuda::managed::array<complex, 2> const a_gpu = a;
//       multi::cuda::managed::array<complex, 2> const b_gpu = b;
//       {
//           multi::cuda::managed::array<complex, 2> c_gpu({2, 2});
//           blas::gemm(ctx, 1., blas::H(a_gpu), b_gpu, 0., c_gpu); // c=ab, c⸆=b⸆a⸆
//           BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//       }
//       {
//           auto c_gpu =+ blas::gemm(&ctx, blas::H(a_gpu), b_gpu);
//           BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//       }
//   }
// }

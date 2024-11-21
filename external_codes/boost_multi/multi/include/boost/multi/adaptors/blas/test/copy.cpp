// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/blas/copy.hpp>  // for copy, copy_n
#include <boost/multi/array.hpp>               // for array, layout_t, subarray

#if defined(NDEBUG)  //  && !defined(__NVCC__) && !(defined(__clang__) && defined(__CUDA__))
	#include <algorithm>  // for transform
	#include <chrono>     // NOLINT(build/c++11) for duration, high_resolution...
	#if __has_include(<execution>) && !defined(__NVCC__) && !defined(__NVCOMPILER)
		#if !((defined(__clang__) && !defined(__apple_build_version__)) && defined(__CUDA__))
			#if (!defined(__INTEL_LLVM_COMPILER) || (__INTEL_LLVM_COMPILER > 20240000))
				#include <execution>  // for execution_policy
			#endif
		#endif
	#endif
	#include <functional>  // for invoke  // IWYU pragma: keep
	#include <iostream>    // for basic_ostream, endl, cout
#endif

#include <complex>   // for operator*, operator+
#include <iterator>  // for size

namespace multi = boost::multi;
namespace blas  = multi::blas;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE)

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(multi_blas_copy_n) {
		multi::array<double, 1> const x = {1.0, 2.0, 3.0, 4.0};  // NOLINT(readability-identifier-length) BLAS naming
		multi::array<double, 1>       y = {5.0, 6.0, 7.0, 8.0};  // NOLINT(readability-identifier-length) BLAS naming
		blas::copy_n(x.begin(), x.size(), y.begin());
		BOOST_TEST( y == x );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_copy) {
		multi::array<double, 1> const x = {1.0, 2.0, 3.0, 4.0};  // NOLINT(readability-identifier-length) BLAS naming
		{
			multi::array<double, 1> y = {5.0, 6.0, 7.0, 8.0};  // NOLINT(readability-identifier-length) BLAS naming
			blas::copy(x, y);                                  // segmentation fault in clang-11
			BOOST_TEST( y == x );
		}
		{
			multi::array<double, 1> y = {5.0, 6.0, 7.0, 8.0};  // NOLINT(readability-identifier-length) BLAS naming
			BOOST_TEST( size(y) == size(x) );
			y() = blas::copy(x);
			BOOST_TEST( y == x );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_real) {
		namespace blas = multi::blas;

		multi::array<double, 2> arr = {
			{1.0,  2.0,  3.0,  4.0},
			{5.0,  6.0,  7.0,  8.0},
			{9.0, 10.0, 11.0, 12.0},
		};

		BOOST_TEST( arr[0][2] ==  3.0 );
		BOOST_TEST( arr[2][2] == 11.0 );

		blas::copy(arr[0], arr[2]);
		BOOST_TEST( arr[0][2] ==  3.0 );
		BOOST_TEST( arr[2][2] ==  3.0 );

		blas::copy(arr[1]({0, size(arr[1])}), arr[2]({0, size(arr[1])}));
		BOOST_TEST( arr[1][3] == 8.0 );
		BOOST_TEST( arr[2][3] == 8.0 );

		multi::array<double, 1> AR3 = blas::copy(arr.rotated()[3]);  // dcopy
		BOOST_TEST( AR3[1] == arr[1][3] );
	}

	BOOST_AUTO_TEST_CASE(multi_blas_copy_row) {
		multi::array<double, 2> const arr = {
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
			{7.0, 8.0, 9.0},
		};
		multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{3}});  // NOLINT(readability-identifier-length) BLAS naming
		blas::copy(arr.rotated()[0], y);
		BOOST_TEST( y == arr.rotated()[0] );
	}

	BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_copy_complex) {
		using complex = std::complex<double>;

		auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> arr = {
			{1.0 + 3.0 * I,  2.0 + 4.0 * I,  3.0 + 5.0 * I,  4.0 + 6.0 * I},
			{5.0 + 0.0 * I,  6.0 + 0.0 * I,  7.0 + 0.0 * I,  8.0 + 0.0 * I},
			{9.0 + 0.0 * I, 10.0 + 0.0 * I, 11.0 + 0.0 * I, 12.0 + 0.0 * I},
		};
		blas::copy(arr[0], arr[2]);
		BOOST_TEST( arr[0][2] == 3.0 + 5.0*I );
	}

#if defined(__INTEL_LLVM_COMPILER)
	std::cout << __INTEL_LLVM_COMPILER << std::endl;
#endif

#if defined(NDEBUG)
	/* transform copy */ {
		multi::array<double, 2> A2D({10000, 10000}, 55.5);
		auto&&                  A2D_block = A2D({1000, 9000}, {1000, 5000});

		multi::array<double, 2> B2D({10000, 10000}, 66.6);
		auto&&                  B2D_block = ~(~B2D({1000, 9000}, {1000, 9000})).strided(2);

		using std::chrono::high_resolution_clock;
		using std::chrono::duration;

		std::cout
			<< "MULTI assignment\n"
			<< std::invoke([&, start_time = high_resolution_clock::now()] {
				   B2D_block = A2D_block;
				   return duration<double>{high_resolution_clock::now() - start_time};
			   }).count()
			<< '\n';

		BOOST_TEST( A2D_block == B2D_block );

		std::cout << "std::transform BLAS\n"
				  << std::invoke([&, start_time = high_resolution_clock::now()] {
						 std::transform(A2D_block.begin(), A2D_block.end(), B2D_block.begin(), [](auto const& row) { return multi::blas::copy(row); });
						 return duration<double>{high_resolution_clock::now() - start_time};
					 }).count()
				  << '\n';

		BOOST_TEST( A2D_block == B2D_block );

	#if __has_include(<execution>) && !defined(__NVCC__) && !defined(__NVCOMPILER)
	#if !((defined(__clang__) ) && defined(__CUDA__)) && (!defined(__INTEL_LLVM_COMPILER) || (__INTEL_LLVM_COMPILER > 20240000))
		#if(__cplusplus >= 202002L)
		#if !defined(__apple_build_version__)
		std::cout << "std::transform par BLAS\n"
				  << std::invoke([&, start_time = high_resolution_clock::now()] {
						 std::transform(std::execution::par, A2D_block.begin(), A2D_block.end(), B2D_block.begin(), [](auto& row) { return multi::blas::copy(row); });
						 return duration<double>{high_resolution_clock::now() - start_time};
					 }).count()
				  << '\n';

		BOOST_TEST( A2D_block == B2D_block );
		#endif
		std::cout << "std::copy par\n"
				  << std::invoke([&, start_time = high_resolution_clock::now()] {
						 std::copy(std::execution::par, A2D_block.begin(), A2D_block.end(), B2D_block.begin());
						 return duration<double>{high_resolution_clock::now() - start_time};
					 }).count()
				  << '\n';

		std::cout << "std::copy par 2\n"
				  << std::invoke([&, start_time = high_resolution_clock::now()] {
						 std::transform(
							 std::execution::par, A2D_block.begin(), A2D_block.end(), B2D_block.begin(), B2D_block.begin(),
							 [](auto const& row_a, auto&& row_b) -> auto&& {
								 std::copy(std::execution::par_unseq, row_a.begin(), row_a.end(), row_b.begin());
								 return std::forward<decltype(row_b)>(row_b);
							 }
						 );
						 return duration<double>{high_resolution_clock::now() - start_time};
					 }).count()
				  << '\n';

		BOOST_TEST( A2D_block == B2D_block );

		std::cout << "std::copy elements par\n"
				  << std::invoke([&, start_time = high_resolution_clock::now()] {
						 std::copy(std::execution::par_unseq, A2D_block.elements().begin(), A2D_block.elements().end(), B2D_block.elements().begin());
						 return duration<double>{high_resolution_clock::now() - start_time};
					 }).count()
				  << '\n';

		BOOST_TEST( A2D_block == B2D_block );
		#endif
	#endif
	#endif

		std::cout << "std::copy\n"
				  << std::invoke([&, start_time = high_resolution_clock::now()] {
						 std::copy(A2D_block.begin(), A2D_block.end(), B2D_block.begin());
						 return duration<double>{high_resolution_clock::now() - start_time};
					 }).count()
				  << '\n';

		BOOST_TEST( A2D_block == B2D_block );

		std::cout << "Multi element assignment\n"
				  << std::invoke([&, start_time = high_resolution_clock::now()] {
						 B2D_block.elements() = A2D_block.elements();
						 return duration<double>{high_resolution_clock::now() - start_time};
					 }).count()
				  << '\n';

		BOOST_TEST( A2D_block == B2D_block );
	}
#endif

	return boost::report_errors();
}

// Copyright 2022-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust universal copy and assignment"
// #include<boost/test/unit_test.hpp>

#include <boost/multi/array.hpp>
//#include<thrust/system/cuda/memory.h>

#include <boost/multi/adaptors/thrust.hpp>

//#include <multi/adaptors/thrust/fix_complex_traits.hpp>


#include <thrust/complex.h>

namespace multi = boost::multi;
using complex = thrust::complex<double>;

template <typename T>
void doNotOptimize(T const& val) {
  asm volatile("" : : "g"(val) : "memory");
}


#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

// #define BOOST_REQUIRE_CLOSE(X, Y, ToL) BOOST_TEST( std::abs( (X) - (Y) ) < (ToL) )

int main() {
BOOST_AUTO_TEST_CASE(thrust_universal_speed_algo) {

	auto const n = 8000;
	{ //cctor
		auto tick = std::chrono::high_resolution_clock::now();
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});
		cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements()*sizeof(complex), 0);

		auto size = A.num_elements()*sizeof(complex)/1e9;
		std::fill_n(raw_pointer_cast(A.data_elements()), A.num_elements(), complex{1.0});

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		double rate = size/time.count();

	}
	{ //cctor
		auto tick = std::chrono::high_resolution_clock::now();
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});

		std::fill_n(raw_pointer_cast(A.data_elements()), A.num_elements(), complex{1.0});

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

		auto size = A.num_elements()*sizeof(complex)/1e9;
		double rate = size/time.count();

		std::cout<<"no  prefetch+cpu_algo rate = "<< rate <<" GB/s\n";
	}
	{ //cctor
		auto tick = std::chrono::high_resolution_clock::now();
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});
		cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements()*sizeof(complex), 0);

		thrust::fill_n(raw_pointer_cast(A.data_elements()), A.num_elements(), complex{1.0});

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

		auto size = A.num_elements()*sizeof(complex)/1e9;
		double rate = size/time.count();

		std::cout<<"dev prefetch+cpu_algo rate = "<< rate <<" GB/s\n";
	}
	{
		auto tick = std::chrono::high_resolution_clock::now();
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});

		thrust::fill_n(A.data_elements(), A.num_elements(), complex{1.0});

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

		auto size = A.num_elements()*sizeof(complex)/1e9;
		double rate = size/time.count();

		std::cout<<"no  prefetch+gpu_algo rate = "<< rate <<" GB/s\n";
	}
	{
		auto tick = std::chrono::high_resolution_clock::now();
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});
		cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements()*sizeof(complex), 0);

		thrust::fill_n(A.data_elements(), A.num_elements(), complex{1.0});

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

		auto size = A.num_elements()*sizeof(complex)/1e9;
		double rate = size/time.count();

		std::cout<<"dev prefetch+gpu_algo rate = "<< rate <<" GB/s\n";
	}
	{
		auto tick = std::chrono::high_resolution_clock::now();
		multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n}, complex{0.0});
		doNotOptimize(A);
		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

		auto size = A.num_elements()*sizeof(complex)/1e9;
		double rate = size/time.count();

		std::cout<<"fill constructor rate = "<< rate <<" GB/s\n";
	}
}

BOOST_AUTO_TEST_CASE(thrust_run) {
	multi::array<long, 1, thrust::cuda::allocator<long>> A(100);
}

return boost::report_errors();

}
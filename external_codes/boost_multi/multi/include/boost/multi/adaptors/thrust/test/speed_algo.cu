// Copyright 2022-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>
#include <boost/multi/detail/what.hpp>

#include <thrust/complex.h>

#include <boost/core/lightweight_test.hpp>

#include <chrono>

namespace multi = boost::multi;
using complex   = thrust::complex<double>;

template<typename T>
void doNotOptimize(T const& val) {
#if defined(_MSC_VER)
	_ReadWriteBarrier();
	(void)val;
#else
	asm volatile("" : : "g"(val) : "memory");
#endif
}

auto universal_memory_supported() -> bool {
	std::cout << "testing for universal memory supported" << std::endl;
	int d;
	cudaGetDevice(&d);
	int is_cma = 0;
	cudaDeviceGetAttribute(&is_cma, cudaDevAttrConcurrentManagedAccess, d);
	if(is_cma) {
		std::cout << "universal memory is supported" << std::endl;
	} else {
		std::cout << "universal memory is NOT supported" << std::endl;
	}
	return (is_cma == 1) ? true : false;
}

int main() {
	// BOOST_AUTO_TEST_CASE(thrust_universal_speed_algo)
	if(universal_memory_supported()) {
		auto const n = 8000;
		{  // cctor
			auto const tick = std::chrono::high_resolution_clock::now();

			multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});

#if(CUDART_VERSION < 13000)
			// cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements() * sizeof(complex), 0);
			cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements() * sizeof(complex), 0, 0);
// #else
// cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements() * sizeof(complex), cudaMemLocation{cudaMemLocationTypeHost, 0}, 0);
#endif

			auto size = A.num_elements() * sizeof(complex) / 1e9;
			std::fill_n(raw_pointer_cast(A.data_elements()), A.num_elements(), complex{1.0});

			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
			double                        rate = size / time.count();
			(void)rate;
		}
		{  // cctor
			auto tick = std::chrono::high_resolution_clock::now();

			multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});

			std::fill_n(raw_pointer_cast(A.data_elements()), A.num_elements(), complex{1.0});

			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

			auto   size = A.num_elements() * sizeof(complex) / 1e9;
			double rate = size / time.count();

			std::cout << "no  prefetch+cpu_algo rate = " << rate << " GB/s\n";
		}
		{  // cctor
			auto tick = std::chrono::high_resolution_clock::now();

			multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});
#if(CUDART_VERSION < 13000)
			// cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements() * sizeof(complex), 0);
			cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements() * sizeof(complex), 0,  cudaStream_t{});
#endif

			thrust::fill_n(raw_pointer_cast(A.data_elements()), A.num_elements(), complex{1.0});

			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

			auto   size = A.num_elements() * sizeof(complex) / 1e9;
			double rate = size / time.count();

			std::cout << "dev prefetch+cpu_algo rate = " << rate << " GB/s\n";
		}
		{
			auto tick = std::chrono::high_resolution_clock::now();

			multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});

			thrust::fill_n(A.data_elements(), A.num_elements(), complex{1.0});

			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

			auto   size = A.num_elements() * sizeof(complex) / 1e9;
			double rate = size / time.count();

			std::cout << "no  prefetch+gpu_algo rate = " << rate << " GB/s\n";
		}
		{
			auto tick = std::chrono::high_resolution_clock::now();

			multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n});
#if(CUDART_VERSION < 13000)
			cudaMemPrefetchAsync(raw_pointer_cast(A.data_elements()), A.num_elements() * sizeof(complex), 0);
#endif
			thrust::fill_n(A.data_elements(), A.num_elements(), complex{1.0});

			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

			auto   size = A.num_elements() * sizeof(complex) / 1e9;
			double rate = size / time.count();

			std::cout << "dev prefetch+gpu_algo rate = " << rate << " GB/s\n";
		}
		{
			auto tick = std::chrono::high_resolution_clock::now();

			multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> A({n, n}, complex{0.0});
			doNotOptimize(A);
			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;

			auto   size = A.num_elements() * sizeof(complex) / 1e9;
			double rate = size / time.count();

			std::cout << "fill constructor rate = " << rate << " GB/s\n";
		}
	}

	// BOOST_AUTO_TEST_CASE(thrust_run)
	{
		multi::array<long, 1, thrust::cuda::allocator<long>> A(100);
	}

	return boost::report_errors();
}

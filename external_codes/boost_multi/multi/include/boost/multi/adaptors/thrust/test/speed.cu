// Copyright 2023-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

#include <thrust/complex.h>

#include <chrono>

namespace multi = boost::multi;

template<>
inline constexpr bool multi::force_element_trivial_default_construction<thrust::complex<double>> = true;

template<>
inline constexpr bool multi::force_element_trivial_default_construction<thrust::complex<float>> = true;

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

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

	// BOOST_AUTO_TEST_CASE(warmup)
	if(universal_memory_supported())
	{
		using T = double;

		auto const n = 8000;

		multi::array<T, 2, thrust::cuda::universal_allocator<T>> src({n, n});
		multi::array<T, 2, thrust::cuda::universal_allocator<T>> dst(extents(src));

		auto const threshold = 0.30;

		auto const size = src.num_elements() * sizeof(T) / 1e9;

		auto const dummy = std::invoke([&] {
			auto start_time = std::chrono::high_resolution_clock::now();
			cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
			auto                          rate = size / time.count();
			// std::cout<<"memcpy    rate = "<< rate <<" GB/s (ratio = 1)\n";
			return rate;
		});
		(void)dummy;

		auto const memcpy_rate = std::invoke([&] {
			auto start_time = std::chrono::high_resolution_clock::now();
			cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
			auto                          rate = size / time.count();
			// std::cout<<"memcpy    rate = "<< rate <<" GB/s (ratio = 1)\n";
			return rate;
		});

		{  // cctor
			auto tick = std::chrono::high_resolution_clock::now();

			auto dst2 = src;

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;

			// std::cout<<"cctor      rate = "<< rate <<" GB/s (ratio = "<< ratio <<")\n";
			if(ratio >= threshold) {
				std::cout << "x";
			}
		}
		{  // assign
			auto tick = std::chrono::high_resolution_clock::now();

			dst = src;

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;

			// std::cout << "assign     rate = "<< rate <<" GB/s (ratio = "<< ratio <<")\n";
			if(ratio >= threshold) {
				std::cout << "x";
			}
		}
		{  // subarray assign
			auto tick = std::chrono::high_resolution_clock::now();

			dst({0, n - 2}, {0, n - 2}) = src({2, n}, {2, n});

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;
			// std::cout << "subasssign rate = "<< rate <<" GB/s (ratio = "<< ratio << ")\n";
			if(ratio >= threshold) {
				std::cout << "x";
			}
		}
	}

	// BOOST_AUTO_TEST_CASE(thrust_nonuniversal_speed)
	{
		using T = ::thrust::complex<double>;
		std::cout << typeid(T).name() << " ******************************************\n";

		auto const n = 8000;

		using AllocatorT = thrust::cuda::allocator<T>;

		multi::array<T, 2, AllocatorT> src({n, n});
		multi::array<T, 2, AllocatorT> dst(extents(src));

		auto const threshold = 0.10;

		auto const size = src.num_elements() * sizeof(T) / 1e9;

		auto const dummy = std::invoke([&] __host__ {
			auto start_time = std::chrono::high_resolution_clock::now();
			cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
			auto                          rate = size / time.count();
			std::cout << "memcpy    rate = " << rate << " GB/s (warmup)\n";
			return rate;
		});
		(void)dummy;

		auto const memcpy_rate = std::invoke([&] __host__ {
			auto start_time = std::chrono::high_resolution_clock::now();
			cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
			auto                          rate = size / time.count();
			std::cout << "memcpy    rate = " << rate << " GB/s (ratio = 1)\n";
			return rate;
		});

		{  // cctor
			auto tick = std::chrono::high_resolution_clock::now();

			auto dst2 = src;

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;

			std::cout << "cctor      rate = " << rate << " GB/s (ratio = " << ratio << ")\n";
			BOOST_TEST(ratio >= threshold);
		}
		{  // assign
			auto tick = std::chrono::high_resolution_clock::now();

			dst = src;

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;

			std::cout << "assign     rate = " << rate << " GB/s (ratio = " << ratio << ")\n";
			BOOST_TEST(ratio >= threshold);
		}
		{  // subarray assign
			auto tick = std::chrono::high_resolution_clock::now();

			dst({0, n - 2}, {0, n - 2}) = src({2, n}, {2, n});

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;

			std::cout << "subasssign rate = " << rate << " GB/s (ratio = " << ratio << ")\n";
			BOOST_TEST(ratio >= threshold);
		}
	}

	// BOOST_AUTO_TEST_CASE(thrust_universal_speed)
	if(universal_memory_supported())
	{
		using T = ::thrust::complex<double>;
		std::cout << typeid(T).name() << " ******************************************\n";

		auto const n = 8000;

		using AllocatorT = thrust::cuda::universal_allocator<T>;

		multi::array<T, 2, AllocatorT> src({n, n});
		multi::array<T, 2, AllocatorT> dst(extents(src));

		auto const threshold = 0.10;

		auto const size = src.num_elements() * sizeof(T) / 1e9;

		auto const dummy = std::invoke([&] __host__ {
			auto start_time = std::chrono::high_resolution_clock::now();
			cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
			auto                          rate = size / time.count();
			std::cout << "memcpy    rate = " << rate << " GB/s (warmup)\n";
			return rate;
		});
		(void)dummy;

		auto const memcpy_rate = std::invoke([&] __host__ {
			auto start_time = std::chrono::high_resolution_clock::now();
			cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), cudaMemcpyDeviceToDevice);
			cudaDeviceSynchronize();
			std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
			auto                          rate = size / time.count();
			std::cout << "memcpy    rate = " << rate << " GB/s (ratio = 1)\n";
			return rate;
		});

		{  // cctor
			auto tick = std::chrono::high_resolution_clock::now();

			auto dst2 = src;

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;

			std::cout << "cctor      rate = " << rate << " GB/s (ratio = " << ratio << ")\n";
			// BOOST_WARN(ratio >= threshold);
		}
		{  // assign
			auto tick = std::chrono::high_resolution_clock::now();

			dst = src;

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;

			std::cout << "assign     rate = " << rate << " GB/s (ratio = " << ratio << ")\n";
			BOOST_TEST(ratio >= threshold);
		}
		{  // subarray assign
			auto tick = std::chrono::high_resolution_clock::now();

			dst({0, n - 2}, {0, n - 2}) = src({2, n}, {2, n});

			std::chrono::duration<double> time  = std::chrono::high_resolution_clock::now() - tick;
			double                        rate  = size / time.count();
			double                        ratio = rate / memcpy_rate;
			std::cout << "subasssign rate = " << rate << " GB/s (ratio = " << ratio << ")\n";
			BOOST_TEST(ratio >= threshold);
		}
	}

	return boost::report_errors();
}

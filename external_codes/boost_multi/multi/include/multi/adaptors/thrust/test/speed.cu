#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust universal copy and assignment"
#include<boost/test/unit_test.hpp>

#include <multi/array.hpp>
#include <multi/adaptors/thrust.hpp>

#include <thrust/complex.h>

namespace multi = boost::multi;
using complex = thrust::complex<double>;

BOOST_AUTO_TEST_CASE(thrust_universal_speed) {

	auto const n = 8000;

	multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> src({n, n});
	multi::array<complex, 2, thrust::cuda::universal_allocator<complex>> dst(extensions(src));

	auto const threshold = 0.2;

	auto const size = src.num_elements()*sizeof(complex)/1e9;

	auto const dummy = std::invoke([&]{
		auto start_time = std::chrono::high_resolution_clock::now();
		cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements()*sizeof(complex), cudaMemcpyDeviceToDevice);
		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
		auto rate = size/time.count();
		std::cout<<"memcpy    rate = "<< rate <<" GB/s (ratio = 1)\n";
		return rate;
	});

	auto const memcpy_rate = std::invoke([&]{
		auto start_time = std::chrono::high_resolution_clock::now();
		cudaMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements()*sizeof(complex), cudaMemcpyDeviceToDevice);
		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
		auto rate = size/time.count();
		std::cout<<"memcpy    rate = "<< rate <<" GB/s (ratio = 1)\n";
		return rate;
	});

	{ //cctor
		auto tick = std::chrono::high_resolution_clock::now();

		auto dst2 = src;

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		double rate = size/time.count();
		double ratio = rate/memcpy_rate;

		std::cout<<"cctor      rate = "<< rate <<" GB/s (ratio = "<< ratio <<")\n";
		BOOST_TEST(ratio >= threshold);
	}
	{ //assign
		auto tick = std::chrono::high_resolution_clock::now();

		dst = src;

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		double rate = size/time.count();
		double ratio = rate/memcpy_rate;

		std::cout << "assign     rate = "<< rate <<" GB/s (ratio = "<< ratio <<")\n";
		BOOST_TEST(ratio >= threshold);
	}
	{ //subarray assign
		auto tick = std::chrono::high_resolution_clock::now();

		dst({0, n - 2}, {0, n - 2}) = src({2, n}, {2, n});

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		double rate = size/time.count();
		double ratio = rate/memcpy_rate;
		std::cout << "subasssign rate = "<< rate <<" GB/s (ratio = "<< ratio << ")\n";
		BOOST_TEST(ratio >= threshold);
	}
}


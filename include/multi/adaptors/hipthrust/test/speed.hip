#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust universal copy and assignment"
#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>

// #include <multi/adaptors/thrust.hpp>

// #include <thrust/complex.h>
//#include <thrust/device_ptr.h>
#include <thrust/system/hip/memory.h>
// #include <thrust/uninitialized_copy.h>

// #include <thrust/universal_allocator.h>

#include <multi/adaptors/thrust/fix_pointer_traits.hpp>

#include <thrust/complex.h>

#include <boost/mpl/list.hpp>

namespace multi = boost::multi;

// template<>
// inline constexpr bool multi::force_element_trivial_default_construction<std::complex<double>> = false;

// template<>
// inline constexpr bool multi::force_element_trivial_default_construction<thrust::complex<double>> = false;

// template<>
// inline constexpr bool multi::force_element_trivial_default_construction<std::complex<float>> = false;

// template<>
// inline constexpr bool multi::force_element_trivial_default_construction<thrust::complex<float>> = false;

using test_types = boost::mpl::list<
	char, unsigned, int,
	::thrust::complex<double>, std::complex<double>,
	::thrust::complex<float>, std::complex<float>,
	double, float>;

BOOST_AUTO_TEST_CASE(warmup) {
	using T = double;

	auto const n = 8000;

	multi::array<T, 2, thrust::hip::universal_allocator<T>> src({n, n});
	multi::array<T, 2, thrust::hip::universal_allocator<T>> dst(extensions(src));

	auto const threshold = 0.30;

	auto const size = src.num_elements() * sizeof(T) / 1e9;

	auto const dummy = std::invoke([&] {
		auto start_time = std::chrono::high_resolution_clock::now();
		auto r = hipMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), hipMemcpyDeviceToDevice);
		assert(r == hipSuccess);
		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
		auto                          rate = size / time.count();
		// std::cout<<"memcpy    rate = "<< rate <<" GB/s (ratio = 1)\n";
		return rate;
	});

	auto const memcpy_rate = std::invoke([&] {
		auto start_time = std::chrono::high_resolution_clock::now();
		auto r = hipMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), hipMemcpyDeviceToDevice);
		assert(r == hipSuccess);
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

BOOST_AUTO_TEST_CASE_TEMPLATE(thrust_universal_speed, T, test_types) {
	std::cout << typeid(T).name() << " ******************************************\n";

	auto const n = 8000;

	multi::array<T, 2, thrust::hip::universal_allocator<T>> src({n, n});
	multi::array<T, 2, thrust::hip::universal_allocator<T>> dst(extensions(src));

	auto const threshold = 0.10;

	auto const size = src.num_elements() * sizeof(T) / 1e9;

	auto const dummy = std::invoke([&] {
		auto start_time = std::chrono::high_resolution_clock::now();
		auto r = hipMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), hipMemcpyDeviceToDevice);
		assert(r == hipSuccess);
		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - start_time;
		auto                          rate = size / time.count();
		std::cout << "memcpy    rate = " << rate << " GB/s (warmup)\n";
		return rate;
	});

	auto const memcpy_rate = std::invoke([&] {
		auto start_time = std::chrono::high_resolution_clock::now();
		auto r = hipMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), hipMemcpyDeviceToDevice);
		assert(r == hipSuccess);
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

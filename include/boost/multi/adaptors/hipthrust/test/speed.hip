// Copyright 2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #include <boost/test/unit_test.hpp>

#include <boost/multi/array.hpp>

#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/system/hip/memory.h>

#include <boost/multi/adaptors/thrust/fix_pointer_traits.hpp>

#include <thrust/complex.h>

#include <boost/mpl/list.hpp>

namespace multi = boost::multi;

// template<> inline constexpr bool multi::force_element_trivial_default_construction<std::complex<double>> = true;
// template<> inline constexpr bool multi::force_element_trivial_default_construction<std::complex<float >> = true;

// template<> inline constexpr bool multi::force_element_trivial_default_construction<thrust::complex<double>> = true;
// template<> inline constexpr bool multi::force_element_trivial_default_construction<thrust::complex<float >> = true;

using test_types = boost::mpl::list<
	double, ::thrust::complex<double>, std::complex<double>,
	float , ::thrust::complex<float> , std::complex<float >,
	// char,
	unsigned, int
>;

template <class Tp>
inline void DoNotOptimize(Tp&& value) noexcept {
	(hipDeviceSynchronize()==hipSuccess)?void():exit(-1);
	asm volatile("" : "+m"(value) : : "memory");
}

constexpr auto measure = [](auto&& F) {
	auto start_time = std::chrono::high_resolution_clock::now();
	F();
	std::chrono::duration<double> const time = std::chrono::high_resolution_clock::now() - start_time;
	DoNotOptimize(F);
	return time.count();
};

BOOST_AUTO_TEST_CASE_TEMPLATE(allocation_speed, T, test_types) {
	std::cout << typeid(T).name() << " ******************************************\n";

	auto const n = 8000L;

	auto const size = n*n * sizeof(T) / 1e9;

	// static_assert( std::is_trivial_v<T> || multi::force_element_trivial_default_construction<T> );

	auto hipmallocfree0 = measure([&]() {
		char* Buffer;
		if(auto r = hipMalloc((void**)&Buffer , n*n * sizeof(T)); r != hipSuccess) {throw std::runtime_error{""};}
		DoNotOptimize(Buffer);
        if(auto r = hipFree(Buffer ); r != hipSuccess) {throw std::runtime_error{""};}
		DoNotOptimize(Buffer);
	});

	auto hipmallocfree1 = measure([&]() {
		char* Buffer;
		if(hipMalloc((void**)&Buffer , n*n * sizeof(T)) != hipSuccess) {throw std::runtime_error{""};}
		DoNotOptimize(Buffer);
        if(auto r = hipFree(Buffer ); r != hipSuccess) {throw std::runtime_error{""};}
		DoNotOptimize(Buffer);
	});

	auto hipmallocfree2 = measure([&]() {
		char* Buffer;
		if(auto r = hipMalloc((void**)&Buffer , n*n * sizeof(T)); r != hipSuccess) {throw std::runtime_error{""};}
		DoNotOptimize(Buffer);
        if(auto r = hipFree(Buffer ); r != hipSuccess) {throw std::runtime_error{""};}
		DoNotOptimize(Buffer);
	});

	auto const hipmallocfree = hipmallocfree2;

	auto const hipallocator = measure([&]() {
		thrust::hip::allocator<T> alloc;
		auto p = alloc.allocate(n*n);
		DoNotOptimize(p);
		alloc.deallocate(p, n*n);
		DoNotOptimize(p);
	});

	BOOST_TEST( (hipallocator / hipmallocfree) > 0.10 );

	auto const hipmultiarray = measure([&]() {
		multi::array<T, 2, thrust::hip::allocator<T>> buffer({n, n});
		DoNotOptimize(buffer);
	});

	BOOST_TEST( (hipmultiarray / hipmallocfree) > 0.10 );

	auto const hipmultiuniversalarray = measure([&]() {
		multi::array<T, 2, thrust::hip::universal_allocator<T>> buffer({n, n});
		DoNotOptimize(buffer);
	});

	BOOST_TEST( (hipmultiuniversalarray / hipmallocfree) > 0.80 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(thrust_speed, T, test_types) {
	std::cout << typeid(T).name() << " ******************************************\n";

	auto const n = 8000;

	multi::array<T, 2, thrust::hip::allocator<T>> src({n, n}, T{12});
	multi::array<T, 2, thrust::hip::allocator<T>> dst(extensions(src), T{33});

	auto const size = src.num_elements() * sizeof(T) / 1e9;

	auto const hipmemcpy0 = measure([&]() {
		if(hipMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), hipMemcpyDeviceToDevice) != hipSuccess) throw std::runtime_error("");
		DoNotOptimize(dst);
	});

	BOOST_REQUIRE( dst == src );

	DoNotOptimize(src);
	DoNotOptimize(dst);

	auto const hipmemcpy1 = measure([&]() {
		if(hipMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), hipMemcpyDeviceToDevice) != hipSuccess) throw std::runtime_error("");
		DoNotOptimize(dst);
	});

	BOOST_REQUIRE( dst == src );

	DoNotOptimize(src);
	DoNotOptimize(dst);

	auto const hipmemcpy = measure([&]() {
		hipMemcpy(raw_pointer_cast(dst.data_elements()), raw_pointer_cast(src.data_elements()), src.num_elements() * sizeof(T), hipMemcpyDeviceToDevice) == hipSuccess?:throw std::runtime_error("");
		DoNotOptimize(dst);
	});

	DoNotOptimize(src);
	DoNotOptimize(dst);

	std::cout << "hipmemcpy         rate = " << size/hipmemcpy << " GB/s\n";

	auto const hip_copy = measure([&]() {
		thrust::copy_n(src.data_elements(), src.num_elements(), dst.data_elements());
		DoNotOptimize(dst);
	});

	DoNotOptimize(src); 
	DoNotOptimize(dst);

	std::cout << "hip_copy          rate = " << size/hip_copy << " GB/s\n";

	BOOST_TEST( hipmemcpy > 0.8* hip_copy );

	auto const multi_copy_assign = measure([&]() {
		dst = src;
		DoNotOptimize(dst);
	});

	DoNotOptimize(src);
	DoNotOptimize(dst);

	std::cout << "multi_copy_assign rate = " << size/multi_copy_assign << " GB/s\n";

	BOOST_TEST( hipmemcpy > 0.8* multi_copy_assign );


	auto const hip_copy_elements = measure([&]() {
		thrust::copy(src.elements().begin(), src.elements().end(), dst.elements().begin());
		DoNotOptimize(dst);
	});

	DoNotOptimize(src); 
	DoNotOptimize(dst);

	std::cout << "hip_copy_elements rate = " << size/hip_copy_elements << " GB/s\n";

	BOOST_TEST( hipmemcpy > 0.8* hip_copy_elements );

	auto const multi_sub_assign = measure([&]() {
		dst({2, n - 2}, {2, n - 2}) = src({2, n - 2}, {2, n - 2});
		DoNotOptimize(dst);
	});

	DoNotOptimize(src);
	DoNotOptimize(dst);

	std::cout << "multi_sub_assign rate = " << size/multi_sub_assign << " GB/s\n";

	BOOST_TEST( hipmemcpy > 0.5* multi_sub_assign );
}

#if 0
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
#endif

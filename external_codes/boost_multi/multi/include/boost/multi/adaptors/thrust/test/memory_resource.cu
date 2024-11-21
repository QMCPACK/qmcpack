// Copyright 2022-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust memory resource"

#include <boost/multi/array.hpp>
#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/system/cuda/memory.h>  // for cuda_pointer

#include <thrust/mr/device_memory_resource.h>
#include <thrust/mr/disjoint_pool.h>   // for thrust::mr::disjoint_unsynchronized_pool_resource
#include <thrust/mr/disjoint_tls_pool.h>  // for thrust::mr::tls_disjoint_pool
#include <thrust/mr/pool.h>  // for thrust::mr::unsynchronized_pool_resource

// #include <boost/mpl/list.hpp>

#include <memory_resource>
#include <numeric>

namespace multi = boost::multi;

template<class MultiArray>
void do_stuff_with_array(typename MultiArray::allocator_type alloc) {
    MultiArray arr1({5, 10}, 99., alloc);

	BOOST_REQUIRE( arr1[3][7] == 99. );

	MultiArray arr2(alloc);

    arr2 = arr1;

    arr1.swap(arr2);

    arr1.clear();
    arr1.reextent({20, 30});
	BOOST_REQUIRE(arr1.num_elements() == 600);
}

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

#if 0
BOOST_AUTO_TEST_CASE(thrust_host_memory_resource) {

    thrust::mr::new_delete_resource memres;

    {
		using Alloc = thrust::mr::allocator<int, thrust::mr::new_delete_resource>;
        Alloc alloc(&memres);

        do_stuff_with_array<multi::array<int, 2, Alloc>>(alloc);
    }

    {
        // virtual calls will be issued - wrapping in a polymorphic wrapper
        thrust::mr::polymorphic_adaptor_resource<void*> adaptor(&memres);
		using Alloc = thrust::mr::polymorphic_allocator<int, void*>;
        Alloc alloc(&adaptor);

        do_stuff_with_array<multi::array<int, 2, Alloc>>(alloc);
    }

	{
		using Pool = thrust::mr::unsynchronized_pool_resource<thrust::mr::new_delete_resource>;
		auto pool = Pool{&memres};
		{
		    using Alloc = thrust::mr::allocator<int, Pool>;
		    auto alloc = Alloc{&pool};

		    do_stuff_with_array<multi::thrust::host_array<int, 2, Alloc>>(alloc);
		}
	}
}

BOOST_AUTO_TEST_CASE(thrust_device_memory_resource) {
    {
        // use the global device_ptr-flavored device memory resource
        using Resource = thrust::device_ptr_memory_resource<thrust::device_memory_resource>;
        thrust::mr::polymorphic_adaptor_resource<thrust::device_ptr<void>> adaptor(
            thrust::mr::get_global_resource<Resource>()
        );
        using Alloc = thrust::mr::polymorphic_allocator<int, thrust::device_ptr<void>>;
        Alloc alloc(&adaptor);

        do_stuff_with_array<multi::array<int, 2, Alloc>>(alloc);

		multi::array<int, 2, Alloc> arr({10, 10}, &adaptor);
    }

    thrust::mr::new_delete_resource memres;

	using Pool = thrust::mr::unsynchronized_pool_resource<thrust::mr::new_delete_resource>;
    Pool pool(&memres);
    {
        typedef thrust::mr::allocator<int, Pool> Alloc;
        Alloc alloc(&pool);

        do_stuff_with_array<multi::thrust::host_array<int, 2, Alloc>>(alloc);
    }

	using DisjointPool = thrust::mr::disjoint_unsynchronized_pool_resource<
        thrust::mr::new_delete_resource,
        thrust::mr::new_delete_resource
    >;

    DisjointPool disjoint_pool(&memres, &memres);
    {
        typedef thrust::mr::allocator<int, DisjointPool> Alloc;
        Alloc alloc(&disjoint_pool);

        do_stuff_with_array<multi::thrust::host_array<int, 2, Alloc>>(alloc);
    }
}

BOOST_AUTO_TEST_CASE(thrust_universal_memory_resource) {
    {
        // use the global device_ptr-flavored device memory resource
        thrust::mr::polymorphic_adaptor_resource<thrust::cuda::universal_pointer<void>> adaptor(
			thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>()
		);
        using Alloc = thrust::mr::polymorphic_allocator<int, thrust::cuda::universal_pointer<void>>;
        Alloc alloc(&adaptor);

        do_stuff_with_array<multi::array<int, 2, Alloc>>(alloc);

		multi::array<int, 2, Alloc> arr({10, 10}, &adaptor);
    }
}

BOOST_AUTO_TEST_CASE(thrust_universal_memory_resource_global_resource_b) {
    // use the global device_ptr-flavored device memory resource
    auto adaptor = thrust::mr::polymorphic_adaptor_resource<thrust::cuda::universal_pointer<void>>{thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>()};
	multi::array<int, 2, thrust::mr::polymorphic_allocator<int, thrust::cuda::universal_pointer<void>>> arr({10, 10}, &adaptor);
}

BOOST_AUTO_TEST_CASE(thrust_universal_memory_resource_global_resource_c) {
    // use the global device_ptr-flavored device memory resource
	auto adaptor = thrust::mr::disjoint_unsynchronized_pool_resource<thrust::system::cuda::universal_memory_resource, thrust::mr::new_delete_resource>
	(
		thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
		thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
    );

	using Alloc = thrust::mr::allocator<int, thrust::mr::disjoint_unsynchronized_pool_resource<
		thrust::system::cuda::universal_memory_resource,
		thrust::mr::new_delete_resource
	>>;

	multi::thrust::mr::array<
		int, 2,
		thrust::mr::disjoint_unsynchronized_pool_resource<
			thrust::system::cuda::universal_memory_resource,
			thrust::mr::new_delete_resource
		>
	> arr({10, 10}, &adaptor);
}

BOOST_AUTO_TEST_CASE(thrust_universal_memory_resource_global_resource_d) {
    // use the global device_ptr-flavored device memory resource
	auto adaptor = thrust::mr::disjoint_unsynchronized_pool_resource<thrust::system::cuda::universal_memory_resource, thrust::mr::new_delete_resource>
	(
		thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
		thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
    );

	multi::thrust::mr::array<
		int, 2,
		thrust::mr::memory_resource<thrust::cuda::universal_pointer<void>>
	> arr({10, 10}, &adaptor);
}

BOOST_AUTO_TEST_CASE(thrust_universal_memory_resource_global_resource_e) {
    // use the global device_ptr-flavored device memory resource
	auto res = thrust::mr::disjoint_unsynchronized_pool_resource(
		thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
		thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
    );

	multi::thrust::pmr::array<int, 2, thrust::cuda::universal_pointer<void>> arr({10, 10}, &res);
}

BOOST_AUTO_TEST_CASE(thrust_universal_memory_resource_global_resource_f) {
    // use the global device_ptr-flavored device memory resource
	auto res = thrust::mr::disjoint_unsynchronized_pool_resource(
		thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
		thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
    );

//  multi::thrust::pmr::array<int, 2, thrust::cuda::universal_pointer<void>> arr({10, 10}, &res);
	multi::thrust::cuda::pmr::universal_array<int, 2> arr({10, 10}, &res);
}

#if 1
template <class Tp>
inline __attribute__((always_inline)) void DoNotOptimize(Tp const& value) {
  asm volatile("" : : "r,m"(value) : "memory");
}

template <class Tp>
inline __attribute__((always_inline)) void DoNotOptimize(Tp& value) {
#if defined(__clang__)
  asm volatile("" : "+r,m"(value) : : "memory");
#else
  asm volatile("" : "+m,r"(value) : : "memory");
#endif
}

BOOST_AUTO_TEST_CASE(thrust_benchmark) {

	auto count = 50;

	{
		auto tick = std::chrono::high_resolution_clock::now();

		for(int64_t i = 0; i != count; ++i) {
			multi::thrust::universal_array<int64_t, 2> arr({1000 - i%10, 1000 + i%10});
			DoNotOptimize(arr);
		}

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		std::cout<< "normal arrays " << time.count() <<std::endl;
	}
	{
		auto tick = std::chrono::high_resolution_clock::now();

		auto res = *thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>();

		for(int64_t i = 0; i != count; ++i) {
			multi::thrust::cuda::pmr::universal_array<int64_t, 2> arr({1000 - i%10, 1000 + i%10}, &res);
			DoNotOptimize(arr);
		}

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		std::cout<< "default resource " << time.count() <<std::endl;
	}
	{
		auto tick = std::chrono::high_resolution_clock::now();

		auto res = thrust::mr::disjoint_unsynchronized_pool_resource(
			thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
			thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
		);

		for(int64_t i = 0; i != count; ++i) {
			multi::thrust::cuda::pmr::universal_array<int64_t, 2> arr({1000 - i%10, 1000 + i%10}, &res);
			DoNotOptimize(arr);
		}

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		std::cout<< "polymorphic pool resource " << time.count() <<std::endl;
	}
	{
		auto tick = std::chrono::high_resolution_clock::now();

		auto res = thrust::mr::disjoint_unsynchronized_pool_resource(
			thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
			thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
		);

		for(int64_t i = 0; i != count; ++i) {
			multi::thrust::mr::array<int, 2, decltype(res)> arr({1000 - i%10, 1000 + i%10}, &res);
			DoNotOptimize(arr);
		}

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		std::cout<< "static pool resource " << time.count() <<std::endl;
	}
	{
		auto tick = std::chrono::high_resolution_clock::now();

		auto res = thrust::mr::disjoint_unsynchronized_pool_resource(
			thrust::mr::get_global_resource<thrust::universal_memory_resource>(),
			thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
		);

		for(int64_t i = 0; i != count; ++i) {
			multi::array<int, 2, thrust::mr::allocator<int, decltype(res)> > arr({1000 - i%10, 1000 + i%10}, &res);
			DoNotOptimize(arr);
		}

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		std::cout<< "2 static pool resource " << time.count() <<std::endl;
	}

}

auto& tls_pool(std::pmr::memory_resource* upstream) {
    static thread_local auto adaptor = std::pmr::unsynchronized_pool_resource(std::pmr::new_delete_resource());
    return adaptor;
}

template<
	class T,
	class Base_ = thrust::mr::allocator<T, thrust::mr::memory_resource<thrust::cuda::universal_pointer<void>>>
//  = std::pmr::polymorphic_allocator<T>
>
struct caching_allocator : Base_ {
	caching_allocator() : Base_{
		&thrust::mr::tls_disjoint_pool(thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(), thrust::mr::get_global_resource<thrust::mr::new_delete_resource>())
	//  &            tls_pool         (std::pmr::new_delete_resource())
	} {}
	caching_allocator(caching_allocator const&) : caching_allocator{} {}
	template<class U> struct rebind {using other = caching_allocator<U>;};
};

BOOST_AUTO_TEST_CASE(thrust_benchmark_contd) {

	auto count = 50;

	{
		auto tick = std::chrono::high_resolution_clock::now();

		for(int64_t i = 0; i != count; ++i) {
			multi::array<int, 2, caching_allocator<int>> arr({1000 - i%10, 1000 + i%10});
		//  auto arr2 = arr;
		//  arr2 = arr;
		//  arr2 = std::move(arr);
		//  DoNotOptimize(arr2);
			DoNotOptimize(arr);
		}

		std::chrono::duration<double> time = std::chrono::high_resolution_clock::now() - tick;
		std::cout<< "caching allocator " << time.count() <<std::endl;
	}
}
#endif

#endif
return boost::report_errors();

}
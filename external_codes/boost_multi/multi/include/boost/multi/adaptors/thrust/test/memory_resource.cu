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

return boost::report_errors();

}
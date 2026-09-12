// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array_ref.hpp>
#include <thrust/system/hip/memory.h>
#include <hip/hip_runtime.h>
#include <cstdio>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <boost/core/lightweight_test.hpp>

auto main() -> int {
	int *r1 = nullptr, *r2 = nullptr;

	using dev_ptr = thrust::pointer<int const, thrust::hip_rocprim::tag, int const&, thrust::use_default>;

	hipMalloc(&r1, 10 * sizeof(int))==hipSuccess?:throw std::runtime_error("hipMalloc r1");;
	hipMalloc(&r2, 10 * sizeof(int))==hipSuccess?:throw std::runtime_error("hipMalloc r2");;

	int host[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

	hipMemcpy(r1, host, sizeof(host), hipMemcpyHostToDevice)==hipSuccess?:throw std::runtime_error("hipMemcpy r1");;
	hipMemcpy(r2, host, sizeof(host), hipMemcpyHostToDevice)==hipSuccess?:throw std::runtime_error("hipMemcpy r2");;

	boost::multi::array_ref<int const, 1, dev_ptr> a(dev_ptr(r1), {10});
	boost::multi::array_ref<int const, 1, dev_ptr> b(dev_ptr(r2), {10});

	BOOST_TEST( a == b );

	return boost::report_errors();
}

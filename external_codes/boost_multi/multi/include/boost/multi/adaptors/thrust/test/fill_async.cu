// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi thrust::fill with async (stream) execution policy"

#include <boost/multi/array.hpp>
#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/fill.h>
#include <thrust/system/cuda/execution_policy.h>  // for thrust::cuda::par_nosync

#include <boost/core/lightweight_test.hpp>

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(bugprone-exception-escape)
#if THRUST_VERSION >= 101600
	cudaStream_t stream;  // NOLINT(cppcoreguidelines-init-variables)
	BOOST_TEST( cudaStreamCreate(&stream) == cudaSuccess );

	multi::thrust::cuda::array<double, 2> darr({1024, 1024});

	// `par_nosync` is the asynchronous execution policy: the algorithm is
	// enqueued on `stream` and the call returns without synchronizing.
	thrust::fill(
		thrust::cuda::par_nosync.on(stream),
		darr.elements().begin(), darr.elements().end(),
		42.0
	);

	// nothing is guaranteed complete until the stream is drained
	BOOST_TEST( cudaStreamSynchronize(stream) == cudaSuccess );

	multi::array<double, 2> host = darr;  // copy back to verify
	BOOST_TEST( host[0][0] == 42.0 );
	BOOST_TEST( host[1023][1023] == 42.0 );

	BOOST_TEST( cudaStreamDestroy(stream) == cudaSuccess );
#endif

	return boost::report_errors();
}

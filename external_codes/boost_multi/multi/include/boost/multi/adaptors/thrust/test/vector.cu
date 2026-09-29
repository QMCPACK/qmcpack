// Copyright 2021-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <boost/core/lightweight_test.hpp>

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

// BOOST_AUTO_TEST_CASE(vector)
{
	// H has storage for 4 integers
	thrust::host_vector<int> H(4);

	// initialize individual elements
	H[0] = 14;
	H[1] = 20;
	H[2] = 38;
	H[3] = 46;

	// H.size() returns the size of vector H
	BOOST_TEST( H.size() == 4 );

	// print contents of H
	BOOST_TEST( H[2] == 38 );

	// resize H
	H.resize(2);

	BOOST_TEST( H.size() == 2 );

	// Copy host_vector H to device_vector D
	thrust::device_vector<int> D = H;

//  f(D.data());

	// elements of D can be modified
	D[0] = 99;
	D[1] = 88;

	// thurst::device_ptr<int> p = D.data();  // doesn't work with CUDA 11.8
	thrust::cuda::pointer<int> p = D.data();  // this works with thrust from CUDA 12.1
	BOOST_TEST( p[0] == 99 );

	BOOST_TEST( D[1] == 88 );
}

return boost::report_errors();

}

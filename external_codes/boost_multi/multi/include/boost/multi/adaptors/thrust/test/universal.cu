// Copyright 2023-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

#include <thrust/system/cuda/memory.h>

// #include <boost/mpl/list.hpp>

#include <iostream>
#include <memory_resource>
#include <numeric>

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
	return (is_cma == 1)?true:false;
}

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

	// BOOST_AUTO_TEST_CASE(thrust_universal_ptr)
	if(universal_memory_supported())
	{
		multi::array<double, 2> Host({1024, 1024});
		std::iota(Host.data_elements(), Host.data_elements() + Host.num_elements(), 12.0);

		multi::array<double, 2, thrust::cuda::universal_allocator<double>> Univ({1024, 1024});

		Univ({0, 10}, {0, 20}) = Host({0, 10}, {0, 20});

		multi::array<double, 2> Hos2({1024, 1024});
		Hos2({0, 10}, {0, 20}) = Univ({0, 10}, {0, 20});

		BOOST_TEST( std::abs( Hos2[0][0] - 12.0 ) < 1e-10 );
	}

	// BOOST_AUTO_TEST_CASE(thrust_universal_ptr_initializer_list)
	if(universal_memory_supported())
	{
		multi::array<double, 1> Host = {1.0, 2.0, 3.0};
		BOOST_TEST( Host.size() == 3 );
		{
			multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ(multi::extents_t<1>{3});
			Univ[0] = 3.0;
			Univ[1] = 2.0;
			Univ[2] = 1.0;

			Host() = Univ();

			BOOST_TEST( Host[0] == 3.0 );
		}
		{
			multi::array<double, 1> tmp = {
				3.0,
				2.0,
				1.0,
			};
			multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ{multi::extents_t<1>{3}};
			Univ = tmp;

			Host() = Univ();

			BOOST_TEST( Host[0] == 3.0 );
		}
		{
			multi::array<double, 1> tmp = {
				3.0,
				2.0,
				1.0,
			};
			multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ{tmp};

			Host() = Univ();

			BOOST_TEST( Host[0] == 3.0 );
		}
		{
			multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ = {
				3.0,
				2.0,
				1.0,
			};

			Host() = Univ();

			std::cout << "host 0 " << Host[0] << '\n';
			BOOST_TEST( Host[0] == 3.0 );
		}
	}

	return boost::report_errors();
}

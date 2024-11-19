// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA thrust universal"
// #include<boost/test/unit_test.hpp>

#include <boost/multi/array.hpp>

#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/system/cuda/memory.h>

#include <boost/mpl/list.hpp>

#include <memory_resource>
#include <numeric>

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

BOOST_AUTO_TEST_CASE(thrust_universal_ptr) {
	multi::array<double, 2                                           > Host({1024, 1024});
	std::iota(Host.data_elements(), Host.data_elements() + Host.num_elements(), 12.0);

	multi::array<double, 2, thrust::cuda::universal_allocator<double>> Univ({1024, 1024});

	Univ({0, 10}, {0, 20}) = Host({0, 10}, {0, 20});

	multi::array<double, 2                                           > Hos2({1024, 1024});
	Hos2({0, 10}, {0, 20}) = Univ({0, 10}, {0, 20});

	BOOST_TEST( Hos2[0][0] == 12.0 );
}

BOOST_AUTO_TEST_CASE(thrust_universal_ptr_initializer_list) {
	multi::array<double, 1                                           > Host = {1.0, 2.0, 3.0};
	BOOST_TEST( Host.size() == 3 );
	{
		multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ(multi::extensions_t<1>{3});
		Univ[0] = 3.0;
		Univ[1] = 2.0;
		Univ[2] = 1.0;

		Host() = Univ();

		BOOST_TEST( Host[0] == 3.0 );
	}
	{
		multi::array<double, 1> tmp = {3.0, 2.0, 1.0,};
		multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ{multi::extensions_t<1>{3}};
		Univ = tmp;

		Host() = Univ();

		BOOST_TEST( Host[0] == 3.0 );
	}
	{
		multi::array<double, 1> tmp = {3.0, 2.0, 1.0,};
		multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ{tmp};

		Host() = Univ();

		BOOST_TEST( Host[0] == 3.0 );
	}
	{
		multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ = {3.0, 2.0, 1.0,};

		Host() = Univ();

		BOOST_TEST( Host[0] == 3.0 );
	}
//  what( thrust::cuda::universal_allocator<double>{} );
//  {
//      multi::array<double, 1, thrust::cuda::universal_allocator<double>> Univ = {3., 2., 1.};

//      Host() = Univ();

//      BOOST_TEST( Host[0] == 3. );
//  }
}

return boost::report_errors();

}
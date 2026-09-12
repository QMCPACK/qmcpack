// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

#include <thrust/complex.h>

#include <boost/core/lightweight_test.hpp>

namespace multi = boost::multi;

class nonbuiltin {
	multi::index val_;
	public:
	nonbuiltin(multi::index val) : val_{val} {}
	__host__ __device__ constexpr nonbuiltin(nonbuiltin const& other) : val_{other.val_} {}  // make it non-trivially copyable
	__host__ __device__ constexpr auto val() const -> multi::index { return val_; }
};

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		multi::array<multi::index, 1> const arr 
			= [] (multi::index i) { return i; } ^ multi::extents_t(10);
		BOOST_TEST( arr[3] == 3 );
    }
    {
		multi::array<multi::index, 1, thrust::cuda::allocator<multi::index>> const arr 
			= [] __host__ __device__ (multi::index i) { return i; } ^ multi::extents_t(10);
		BOOST_TEST( arr[3] == 3 );
	}
    // {
	// 	multi::array<multi::index, 1, thrust::cuda::allocator<multi::index>> const arr 
	// 		= [] __device__ (multi::index i) { return i; } ^ multi::extents_t(10);
	// 	BOOST_TEST( arr[3] == 3 );
	// }
    {
		multi::array<multi::index, 2> const arr 
			= [] (multi::index i, multi::index j) { return i + j; } ^ multi::extents_t(10, 10);
		BOOST_TEST( arr[3][4] == 3 + 4 );
	}
	{
		multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr 
			= [] __host__ __device__(multi::index i, multi::index j) { return i + j; } ^ multi::extents_t(10, 10);
		BOOST_TEST( arr[3][4] == 3 + 4 );
	}
	{
		multi::index a = 99;
		multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr 
			= [a] __host__ __device__(multi::index i, multi::index j) { return i + j + a; } ^ multi::extents_t(10, 10);
		BOOST_TEST( arr[3][4] == 3 + 4 + 99);
	}
	{
		nonbuiltin nbi(99);
		multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr 
			= [nbi] __host__ __device__(multi::index i, multi::index j) { return i + j + nbi.val(); } ^ multi::extents_t(10, 10);
		BOOST_TEST( arr[3][4] == 3 + 4 + 99);
	}
	// {
	// 	nonbuiltin nbi(99);
	// 	multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr 
	// 		= [nbi] __host__ __device__(auto i, auto j) { return i + j + nbi.val(); } ^ multi::extents_t(10, 10);
	// 	BOOST_TEST( arr[3][4] == 3 + 4 + 99);
	// }
	{
		nonbuiltin nbi(99);
		multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr 
			= [nbi] (multi::index i, multi::index j) { return i + j + nbi.val(); } ^ multi::extents_t(10, 10);
		BOOST_TEST( arr[3][4] == 3 + 4 + 99);
	}
	// {
	// 	auto fun = [] __device__(int i, int j) { return i + j; };
	// 	multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr 
	// 		= device{fun} ^ multi::extents_t(10, 10);
	// 	BOOST_TEST( arr[3][4] == 3 + 4 + 99);
	// }
	{
		multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr
			= [](multi::index i, multi::index j) constexpr { return i + j; } ^ multi::extents_t(10, 10);
		BOOST_TEST( arr[3][4] == 3 + 4 );
	}
	// {
	// 	multi::array<multi::index, 2, thrust::cuda::allocator<multi::index>> const arr = [](multi::index i, multi::index j) { return i + j; } ^ multi::extents_t(10, 10);
	// 	BOOST_TEST( arr[3][4] == 3 + 4 );
	// }

	// multi::thrust::universal_array<thrust::complex<double>, 2> arr_gold({50, 70});

	// auto [is, js] = arr_gold.extensions();
	// for(auto i : is) {
	// 	for(auto j : js) {
	// 		arr_gold[i][j] = thrust::complex<double>{static_cast<double>(i), static_cast<double>(j)};
	// 	}
	// }

	// multi::thrust::universal_array<thrust::complex<double>, 2> arr =
	// 	[] __host__ __device__ (multi::index i, multi::index j) {
	// 		return thrust::complex<double>{static_cast<double>(i), static_cast<double>(j)};
	// 	} ^
	// 	multi::extents_t<2>({50, 70});

	// BOOST_TEST( arr == arr_gold );

	return boost::report_errors();
}

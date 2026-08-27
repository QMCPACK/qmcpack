// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/generate.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/transform.h>

#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <iostream>
// #include <nvfunctional>

__device__ int square_device(int x) {
	return x * x;
}

struct square_functor {
	__device__ int operator()(int x) const {
		return square_device(x);
	}
};

namespace multi = boost::multi;

// #define _HOST

#ifndef _HOST
template<class T> using vector = multi::thrust::device_array<T, 1>;
#define DEV __device__
#else
template<class T> using vector = multi::thrust::host_array<T, 1>;
#define DEV
#endif

int main() {
	int const N = 8;

	{
		vector<int> d_out(N);

		auto first = thrust::counting_iterator<int>(0);
		auto last  = first + N;

		thrust::transform(
			thrust::device,
			first,
			last,
			d_out.begin(),
			[] DEV(int x) { return x * x; }
		);
	}

	// auto c2 = multi::thrust::device_restriction(multi::extents_t<1>(N), [a = 5] __device__(int x) { return x * x + a; });
	// auto c2 = multi::thrust::device_restriction(multi::extents_t<1>(N), multi::thrust::device_function<int>([a = 5] __device__(int x) { return x * x + a; }));

	// [[maybe_unused]] auto fff = [a = 5] __host__(int x) { return x * x + a; };
	// auto c2
	//     = multi::thrust::device_restriction(multi::extents_t<1>(N),
	//         multi::thrust::device_function<
	//             decltype((fff)(multi::index{}))
	//         >([a = 5] __device__(int x) { return x * x + a; })
	//     );

	// [[maybe_unused]] auto fff = [a = 5] __host__(int x) { return x * x + a; };
	// using result_type = decltype((fff)(multi::index{}));
	// auto c2
	//     = multi::thrust::device_restriction(multi::extents_t<1>(N),
	//         multi::thrust::device_function<
	//             result_type
	//         >([a = 5] __device__(int x) { return x * x + a; })
	//     );

	// auto c2
	//     = multi::thrust::device_restriction(
	//         multi::extents_t<1>(N),
	//         [=]() {
	//             [[maybe_unused]] auto fff = [a = 5] __host__(int x) { return x * x + a; };
	//             using result_type = decltype((fff)(multi::index{}));
	//             return multi::thrust::device_function<result_type>([a = 5] __device__(int x) { return x * x + a; });
	//         }()
	// );

	// auto c2
	// = multi::thrust::device_restriction(
	//         multi::extents_t<1>(N),
	//         [a = 5]() {
	//             // [[maybe_unused]] auto fff = [=] __host__(int x) { return x * x + a; };
	//             // using result_type = typename multi::thrust::result_helper<decltype(fff)>::type;
	//             return [=] __device__(int x) { return x * x + a; };
	//         }()
	// );

// 	[[maybe_unused]] auto fdev = [] __device__(int x) { return x * x + 1.9; };
// 	using rt1                  = std::invoke_result<decltype(fdev), int>::type;
// 	using rt2                  = typename multi::thrust::result_helper<decltype(fdev)>::type;
// 	// multi::detail::what<rt>();
// 	static_assert(std::is_same_v<rt1, double>);
// 	static_assert(std::is_same_v<rt2, double>);

// #if defined(__CUDACC__) &&        \
// 	(__CUDACC_VER_MAJOR__ < 12 || \
// 	 (__CUDACC_VER_MAJOR__ == 12 && __CUDACC_VER_MINOR__ <= 5))

// 	auto c2 = multi::thrust::device_restriction(
// 		multi::extents_t<1>(N),
// 		BOOST_MULTI_DEVICE_LAMBDA_LEGACY(BOOST_MULTI_CAPTURE(a = 5, b = 0), (auto x) { int c = a, d = b; return  x * x + c + d; })
// 	);

// #else

// 	auto c2 = multi::thrust::device_restriction(
// 		multi::extents_t<1>(N),
// 		[ a = 5, b = 0 ] BOOST_MULTI_DEVICE_LAMBDA((auto x) { int c = a, d = b; return  x * x + c + d; })
// 	);

// #endif

// 	auto it = c2.begin();

// 	thrust::copy(c2.begin(), c2.end(), d_out.begin());

// 	multi::thrust::host_array<int, 1> h_out = d_out;

// 	BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
// 	BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );

	// CPU memory and execution, iterator holds function by POINTER semantics
	{
		auto restr = multi::restricted<1>(
			[a = 5, b = 0](auto x) -> int { int c = a, d = b; return static_cast<int>(x * x + c + d); },
			{N}
		);

		multi::thrust::host_array<int, 1> h_out({N}, int{});

		thrust::copy(restr.begin(), restr.end(), h_out.begin());

		BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
		BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );
	}

	// CPU memory and execution, iterator holds function by POINTER semantics
	{
		auto fun   = [a = 5, b = 0](auto x) -> int { int c = a, d = b; return static_cast<int>(x * x + c + d); };
		auto restr = multi::restricted<1>(
			multi::val(fun),
			{N}
		);

		multi::thrust::host_array<int, 1> h_out({N}, int{});

		thrust::copy(restr.begin(), restr.end(), h_out.begin());

		BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
		BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );
	}

	// GPU memory and execution, iterator holds function by VALUE semantics
	{
		auto dev_restr = multi::thrust::device_restricted<1>(
			[a = 5, b = 0] __device__(auto x) { int c = a, d = b; return  x * x + c + d; },
			{N}
		);

		multi::thrust::device_array<int, 1> d_out({N}, int{});

		thrust::copy(
			thrust::device,  // ensure this is run in the GPU (optional)
			dev_restr.begin(), dev_restr.end(), d_out.begin()
		);

		multi::thrust::host_array<int, 1> h_out = d_out;

		BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
		BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );
	}

	// GPU memory and execution, iterator holds function by POINTER semantics
	// THIS DOESN'T WORK AND SHOULDN'T WORK
	// {
	// 	auto device_restr = multi::restricted<1>(
	// 		[a = 5, b = 0] __device__(auto x) { int c = a, d = b; return  x * x + c + d; },
	// 		{N}
	// 	);

	// 	multi::thrust::device_array<int, 1> d_out({N}, int{});

	// 	// vvv this copy fails to actually exectute the device code
	// 	thrust::copy(
	// 		thrust::device,  // ensure this is run in the GPU (optional)
	// 		device_restr.begin(), device_restr.end(), d_out.begin()
	// 	);

	// 	multi::thrust::host_array<int, 1> h_out = d_out;

	// 	// vvv test fails
	// 	BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
	// 	BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );
	// }

	{
		auto dev_restr = multi::thrust::device_restricted<1>(
			multi::val([a = 5, b = 0] __device__(int x) { int c = a, d = b; return  x * x + c + d; }),
			{N}
		);

		multi::thrust::device_array<int, 1> d_out({N}, int{});

		thrust::copy(
			thrust::device,  // ensure this is run in the GPU (optional)
			dev_restr.begin(), dev_restr.end(), d_out.begin()
		);

		multi::thrust::host_array<int, 1> h_out = d_out;

		BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
		BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );
	}
    {
		auto val  = multi::val([a = 5, b = 0](auto x) { int c = a, d = b; return  x * x + c + d; });
		auto pval = &val;
		BOOST_TEST( (*pval)(2) == 2*2 + 5 + 0 );
	}
	// {
	// 	auto device_restr = multi::restricted<1>(
	// 		multi::val([a = 5, b = 0] __device__(int x) { int c = a, d = b; return  x * x + c + d; }),
	// 		{N}
	// 	);

	// 	multi::thrust::device_array<int, 1> d_out({N}, int{999});

	// 	thrust::copy(
	// 		thrust::device,  // ensure this is run in the GPU (optional)
	// 		device_restr.begin(), device_restr.end(), d_out.begin()
	// 	);

	// 	multi::thrust::host_array<int, 1> h_out = d_out;

	// 	std::cout << " h_out[1] = " << h_out[1] << std::endl;
	// 	std::cout << " h_out[2] = " << h_out[2] << std::endl;

	// 	BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
	// 	BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );
	// }

	// GPU memory and execution, iterator holds function by POINTER semantics
	// {
	// 	auto device_restr = multi::restricted<1>(
	// 		[a = 5, b = 0] __device__(auto x) { int c = a, d = b; return  x * x + c + d; },
	// 		{N}
	// 	);

	// 	multi::thrust::device_array<int, 1> d_out({N}, int{});

	// 	// vvv this copy fails to actually exectute the device code
	// 	thrust::copy(
	// 		thrust::device,  // ensure this is run in the GPU (optional)
	// 		device_restr.begin(), device_restr.end(), d_out.begin()
	// 	);

	// 	multi::thrust::host_array<int, 1> h_out = d_out;

	// 	// vvv test fails
	// 	BOOST_TEST( h_out[1] == 1*1 + 5 + 0 );
	// 	BOOST_TEST( h_out[2] == 2*2 + 5 + 0 );
	// }

	{
		thrust::host_vector<int> vec(10);
		std::generate_n(vec.begin(), vec.size(), [a = 1, b = 2]() { return a + b; });

		BOOST_TEST( vec[3] == 1 + 2 );
	}
	{
		::thrust::device_vector<int> vec(10);
		::thrust::generate_n(
			thrust::device,
			vec.begin(), vec.size(),
			[a = 1, b = 2] __device__() { return a + b; }
		);

		thrust::host_vector<int> h_vec = vec;
		;
		BOOST_TEST( h_vec[3] == 1 + 2 );
	}
	// compile error static vars no defined in device code
	// {
	//     ::thrust::device_vector<int> vec(10);
	//     static int a = 1, b = 2;
	//     ::thrust::generate_n(
	//         thrust::device,
	//         vec.begin(), vec.size(),
	//         [] __device__ () {return a + b;}
	//     );

	//     thrust::host_vector<int> h_vec = vec;;
	//     BOOST_TEST( h_vec[3] == 1 + 2 );
	// }
	// {
	//     ::thrust::device_vector<int> vec(10);
	//     int a = 1;
	//     ::thrust::generate_n(
	//         thrust::device,
	//         vec.begin(), vec.size(),
	//         [ap = &a, b = 2] __device__ () {return *ap + b;}
	//     );

	//     thrust::host_vector<int> h_vec = vec;;
	//     BOOST_TEST( h_vec[3] == 1 + 2 );
	// }

	// gives error:   what():  parallel_for: failed to synchronize: cudaErrorIllegalAddress: an illegal memory access was encountered
	// {
	//     thrust::device_vector<int> vec(10);
	//     auto gen = [a = 1, b = 2] __device__ () {return a + b;};
	//     thrust::generate_n(
	//         thrust::device,
	//         vec.begin(), vec.size(),
	//         std::ref(gen)
	//     );

	//     thrust::host_vector<int> h_vec = vec;;
	//     BOOST_TEST( h_vec[3] == 1 + 2 );
	// }
	// gives an error: cudaErrorInvalidPc: invalid program counter
	// {
	//     thrust::device_vector<int> vec(10);
	//     auto gen = [a = 1, b = 2] __device__ () {return a + b;};
	//     thrust::generate_n(
	//         thrust::device,
	//         vec.begin(), vec.size(),
	//         nvstd::function<int()>(gen)
	//     );

	//     BOOST_TEST( vec[3] == 1 + 2 );
	// }

	return boost::report_errors();
}

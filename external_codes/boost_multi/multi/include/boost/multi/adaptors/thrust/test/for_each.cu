// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>
#include <boost/multi/detail/what.hpp>

#include <boost/core/lightweight_test.hpp>

// GCC 12 + nvcc < 12.1 incompatibility: avx512fp16intrin.h uses _Float16 which older nvcc doesn't support
#if defined(TBB_FOUND) && !defined(__NVCC__)
#if !defined(__clang__)
#include <execution>
#include <numeric>
#endif
#endif
// #include <boost/timer/timer.hpp>

#include <chrono>

namespace multi = boost::multi;

namespace {

class auto_timer : std::chrono::high_resolution_clock {
	std::string label_;
	time_point  start_;

 public:
	explicit auto_timer(char const* label) : label_{label} {
		cudaDeviceSynchronize() == cudaSuccess
			? void()
			: assert(0);  // NOLINT(misc-include-cleaner) the header is included
						  // conditionally
		start_ = now();
	}

	auto_timer(auto_timer const&) = delete;
	auto_timer(auto_timer&&)      = delete;

	auto operator=(auto_timer const&) -> auto_timer& = delete;
	auto operator=(auto_timer&&) -> auto_timer&      = delete;

	auto_timer() : auto_timer("") {}

	~auto_timer() {
		cudaDeviceSynchronize() == cudaSuccess ? void() : assert(0);
		auto const count = std::chrono::duration<double>(now() - start_).count();
		std::cerr << label_ << ": " << count << " sec\n";
	}
};
}  // namespace

auto main()
	-> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

	using T = double;
	{
		auto cpu_own = multi::array<T, 3>({64, 64, 64}, 0);

		auto&& cpu = cpu_own();
		{
			auto_timer const _{"std::for_each"};
			std::for_each(cpu.begin(), cpu.end(), [](auto&& plane) {
				for(auto&& row : plane) {     // NOLINT(altera-unroll-loops)
					for(auto&& elem : row) {  // NOLINT(altera-unroll-loops)
						elem += std::sqrt(std::pow(elem, 1.5) + std::sin(elem));
					}
				}
			});
		}

		auto gpu_par = multi::thrust::device_array<T, 3>({64, 64, 64}, 0.0);
		{
			auto_timer t{"nested thrust::for_each(thrust::cuda::par)"};
			for(auto&& plane : gpu_par) {     // NOLINT(altera-unroll-loops)
				for(auto&& row : std::forward<decltype(plane)>(plane)) {     // NOLINT(altera-unroll-loops)
					thrust::for_each(
						thrust::cuda::par,
						row.begin(), row.end(), [] __device__ (auto& elem) {
							elem += std::sqrt(std::pow(elem, 1.5) + std::sin(elem));
						});
				}
			}
		}
		// {
		// 	auto_timer t{"thrust::for_each(thrust::cuda::par)"};
		// 	thrust::for_each(
		// 		thrust::cuda::par,
		// 		gpu_par.begin(), gpu_par.end(), [] __device__ (auto&& plane) {
		// 		for(auto&& row : plane) {     // NOLINT(altera-unroll-loops)
		// 			for(auto&& elem : row) {  // NOLINT(altera-unroll-loops)
		// 				elem += std::sqrt(std::pow(elem, 1.5) + std::sin(elem));
		// 			}
		// 		}
		// 	});
		// }
		{
			auto_timer t{"counting thrust::for_each(thrust::cuda::par)"};
			thrust::for_each(
				thrust::cuda::par,
				thrust::counting_iterator{gpu_par.begin()},
				thrust::counting_iterator{gpu_par.end()},  [] __device__ (auto it_plane) { auto&& plane = *it_plane;
				for(auto&& row : plane) {     // NOLINT(altera-unroll-loops)
					for(auto&& elem : row) {  // NOLINT(altera-unroll-loops)
						elem += std::sqrt(std::pow(elem, 1.5) + std::sin(elem));
					}
				}
			});
		}
		{
			auto_timer t{"iterator thrust::for_each(thrust::cuda::par)"};
			thrust::for_each(thrust::cuda::par, gpu_par.extent().begin(), gpu_par.extent().end(), [gpu_par = gpu_par.begin()] __device__ (auto iplane) {
				auto&& plane = gpu_par[iplane]; 
				for(auto&& row : plane) {     // NOLINT(altera-unroll-loops)
					for(auto&& elem : row) {  // NOLINT(altera-unroll-loops)
						elem += std::sqrt(std::pow(elem, 1.5) + std::sin(elem));
					}
				}
			});
		}
		{
			auto_timer const _{"thrust::for_each(thrust::cuda::par, elements)"};
			thrust::for_each(thrust::cuda::par, gpu_par.elements().begin(), gpu_par.elements().end(), [] __device__(auto& elem) {
				elem += std::sqrt(std::pow(elem, 1.5) + std::sin(elem));
			});
		}
	}

	return boost::report_errors();
}

// Copyright 2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

// #include <thrust/device_allocator.h>
#include <thrust/execution_policy.h>  // Include for execution policies
#include <thrust/transform_reduce.h>  // for thrust::transform_reduce

#include <boost/core/lightweight_test.hpp>

// #include <execution>  // for std::execution::par, doesn't work on gcc 13.3 and nvcc 12.0
#include <iostream>
#include <numeric>  // for std::transform_reduce

namespace multi = boost::multi;

struct v2d {
	double x;
	double y;

	friend __host__ __device__ constexpr auto x(v2d const& self) { return self.x; }
	friend __host__ __device__ constexpr auto y(v2d const& self) { return self.y; }
};

__host__ __device__ constexpr auto v(double dist) -> double {
	return dist * dist;
}

// declare universal array type
// pros: can be used in CPU, and GPU
// cons: members cannot be used directly on array elements (ref wrapped), careless use might cause a lot of page faults (CPU<->GPU page swapping)
template<class T, multi::dimensionality_type D>
using array = multi::thrust::universal_array<T, D>;

// pros: familiarity
// cons: bug prone, lots of state and variables, for init is ugly
template<class Array1D, class Array2D>
auto energy_raw_loops(Array1D const& positions, Array2D const& neighbors) {
	double ret = 0;
	for(multi::index i = 0; i != positions.size(); ++i) {
		auto const& neighbors_i = neighbors[i];
		for(multi::index j = 0; j != neighbors_i.size(); ++j) {
			if(neighbors_i[j] == -1) {
				continue;
			}  // or add zero
			auto dist = std::abs(x(positions[i]) - x(positions[neighbors_i[j]]));
			ret += v(dist);
		}
	}
	return ret;
}

// pros: familiar, more compact
// cons: bug prone
template<class Array1D, class Array2D>
auto energy_range_loops(Array1D const& positions, Array2D const& neighbors) {
	double ret = 0;
	for(auto const i : positions.extent()) {
		auto const positions_i = positions[i];
		for(auto const& nbidx : neighbors[i]) {
			if(nbidx == -1) {
				continue;
			}
			auto dist = std::abs(x(positions_i) - x(positions[nbidx]));
			ret += v(dist);
		}
	}
	return ret;
}

// pros: familiar, less state
// const: too nested
template<class Array1D, class Array2D>
auto energy_reduce_in_loop(Array1D const& positions, Array2D const& neighbors) {
	double ret = 0.0;
	for(auto const i : positions.extent()) {
		ret += std::transform_reduce(
			neighbors[i].begin(), neighbors[i].end(), 0.0, std::plus<>{},
			[positions_i = positions[i], &positions](auto nbidx) {
				return nbidx == -1 ? 0.0 : v(std::abs(x(positions_i) - x(positions[nbidx])));
			}
		);
	}
	return ret;
}

// pros: single ouput, no state, ready for some parallelism, ready for CUDA thrust
// cons: unfamiliar, too many nested structures, return type buried
template<class Array1D, class Array2D>
auto energy_nested_reduce(Array1D const& positions, Array2D const& neighbors) {
	return std::transform_reduce(
		positions.extent().begin(), positions.extent().end(),
		0.0, std::plus<>{},
		[&positions, &neighbors](auto i) {
			return std::transform_reduce(
				neighbors[i].begin(), neighbors[i].end(),
				0.0, std::plus<>{},
				[positions_i = positions[i], &positions](auto nbidx) {
					return nbidx == -1 ? 0.0 : v(std::abs(x(positions_i) - x(positions[nbidx])));
				}
			);
		}
	);
}

// // pros: single ouput, no state, parallel
// // cons: unfamiliar, too many nested structures, parallelization is partial and nested
// auto energy_nested_par_reduce(auto const& positions, auto const& neighbors) {
// 	return std::transform_reduce(
// 		std::execution::par,
// 		positions.extension().begin(), positions.extension().end(),
// 		0.0, std::plus<>{},
// 		[&positions, &neighbors](auto i) {
// 			return std::transform_reduce(
// 				std::execution::unseq,
// 				neighbors[i].begin(), neighbors[i].end(),
// 				0.0, std::plus<>{},
// 				[positions_i = positions[i], &positions](auto nbidx) {
// 					return nbidx == -1 ? 0.0 : v(std::abs(x(positions_i) - x(positions[nbidx])));
// 				}
// 			);
// 		}
// 	);
// }

// pros: correct parallelization, no state
// cons: unfamiliar, needs coordinate decomposition, use special Multi features
// auto energy_flatten_par_reduce(auto const& positions, auto const& neighbors) {
// 	return std::transform_reduce(
// 		std::execution::par,
// 		neighbors.extents().elements().begin(), neighbors.extents().elements().end(),
// 		0.0,
// 		std::plus<>{},
// 		[&positions, &neighbors](multi::array<int, 2>::indexes c) {
// 			auto [i, j] = c;
// 			return neighbors[i][j] == -1 ? 0.0 : v(std::abs(x(positions[i]) - x(positions[neighbors[i][j]])));
// 		}
// 	);
// }

#ifdef _MSC_VER
template<class Positions, class Neighbors>
struct inner_coordinates {
	Positions positions;
	Neighbors neighbors;

	__host__ __device__ inner_coordinates(Positions positions, Neighbors neighbors) : positions{positions}, neighbors{neighbors} {}

	__host__ __device__ auto operator()(array<int, 2>::indexes ij) const -> double {
		auto [i, j] = ij;
		return neighbors[i][j] == -1 ? 0.0 : v(std::abs(x(positions[i]) - x(positions[neighbors[i][j]])));
	}
};
#endif

// pros: GPU optimized, runs completely on GPU
// const: verbose, enclosing function cannot deduce return, or auto parameters, lamba captures need special care, order of arguments is different from STL, may need a complete diffeent algorithm to extract different information: e.g. thrust::reduce_by_key
template<class Arr1D, class Arr2D>
auto energy_flatten_gpu_reduce(Arr1D const& positions, Arr2D const& neighbors) -> double {
	return thrust::transform_reduce(
		thrust::cuda::par,
		neighbors.extents().elements().begin(), neighbors.extents().elements().end(),
#ifdef _MSC_VER
		inner_coordinates{positions.begin(), neighbors.begin()},
#else
		[positions = positions.begin(), neighbors = neighbors.begin()] __device__(array<int, 2>::indexes ij) -> double {
			auto [i, j] = ij;
			return neighbors[i][j] == -1 ? 0.0 : v(std::abs(x(positions[i]) - x(positions[neighbors[i][j]])));
		},
#endif
		0.0,
		std::plus<>{}
	);
}

template<class V2D, class Positions>
struct inner {
	V2D       positions_i;
	Positions positions;

	__host__ __device__ inner(V2D positions_i, Positions positions) : positions_i{positions_i}, positions{positions} {}

	__host__ __device__ auto operator()(array<int, 2>::index nbidx) const {
		return nbidx == -1 ? 0.0
						   : v(std::abs(x(positions_i) - x(positions[nbidx])));
	}
};


template<class Positions, class Neighbors>
struct outer {
	Positions positions;
	Neighbors neighbors;

	__host__ __device__ outer(Positions positions, Neighbors neighbors) : positions{positions}, neighbors{neighbors} {}

	__host__ __device__ auto operator()(array<int, 2>::index i) -> double {
		return thrust::transform_reduce(
			thrust::device,
			neighbors[i].begin(), neighbors[i].end(),
			inner{positions[i], positions},
			0.0,
			std::plus<>{}
		);
	}
};

#if !defined(__clang_major__)
// no pros: this for testing purposed only
// const: requires auxiliary class, and an extra auxiliary class (outer) on Windows
template<class Arr1D, class Arr2D>
auto energy_gpu_nested_reduce(Arr1D const& positions, Arr2D const& neighbors) -> double {
	return thrust::transform_reduce(
		thrust::cuda::par,
		positions.extent().begin(), positions.extent().end(),
	#if !defined(_MSC_VER)
		[positions = positions.begin(), neighbors = neighbors.begin()] __device__(typename Arr1D::index i) -> double {
			return thrust::transform_reduce(
				thrust::device,
				neighbors[i].begin(), neighbors[i].end(),
				inner{positions[i], positions},
				0.0,
				std::plus<>{}
			);
		},
	#else
		outer{positions.begin(), neighbors.begin()},
	#endif
		0.0,
		std::plus<>{}
	);
}

#endif

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
	return (is_cma == 1) ? true : false;
}

auto main() -> int {
	if(universal_memory_supported()) {
		array<v2d, 1> positions = {
			{1.0, 0.0},
			{2.0, 0.0},
			{3.0, 0.0},
			{4.0, 0.0},
			{5.0, 0.0}
		};

		array<array<v2d, 1>::index, 2> neighbors = {
			{1, 2, -1, -1}, /* of at 0*/
			{0, 2,  3, -1}, /* of at 1*/
			{0, 1,  3,  4}, /* of at 2*/
			{1, 2,  4, -1}, /* of at 3*/
			{2, 3, -1, -1}  /* of at 4*/
		};

		{
			auto en = energy_raw_loops(positions, neighbors);
			BOOST_TEST( en == 32.0 );
		}
		{
			auto en = energy_range_loops(positions, neighbors);
			BOOST_TEST( en == 32.0 );
		}
		{
			auto en = energy_reduce_in_loop(positions, neighbors);
			BOOST_TEST( en == 32.0 );
		}
		{
			auto en = energy_nested_reduce(positions, neighbors);
			BOOST_TEST( en == 32.0 );
		}
// {
// 	auto en = energy_nested_par_reduce(positions, neighbors);
// 	BOOST_TEST( en == 32.0 );
// }
// {
// 	auto en = energy_flatten_par_reduce(positions, neighbors);
// 	BOOST_TEST( en == 32.0 );
// }
#if !defined(__clang_major__)
		{
			auto en = energy_flatten_gpu_reduce(positions, neighbors);
			BOOST_TEST( en == 32.0 );
		}
#endif
#if !defined(__clang_major__)
		{  // this is not recommended, it is for testing purposes
			auto en = energy_gpu_nested_reduce(positions, neighbors);
			BOOST_TEST( en == 32.0 );
		}
#endif
	}

	return boost::report_errors();
}

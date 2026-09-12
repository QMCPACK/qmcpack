#define ENABLE_GPU 1

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>

#include <thrust/universal_allocator.h>

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

#include <thrust/reduce.h>
#include <thrust/iterator/discard_iterator.h>

#include <boost/multi/adaptors/thrust/reduce_by_index.hpp>

#include <chrono>

template<class Tp>
inline
#if defined(_MSC_VER)
	__forceinline
#else
	__attribute__((always_inline))
#endif
	void
	DoNotOptimize(Tp const& value) {  // NOLINT(readability-identifier-naming)
#if defined(_MSC_VER)
	_ReadWriteBarrier(); (void)value;
#else
	asm volatile("" : : "r,m"(value) : "memory");  // NOLINT(hicpp-no-assembler)
#endif
}

namespace multi = boost::multi;

template<class kernel_type, class array_type>
__global__ void reduce_kernel_vr(multi::ssize_t sizex, multi::ssize_t sizey, kernel_type kernel, array_type odata) {

	extern __shared__ char shared_mem[];
	auto                   reduction_buffer = (typename array_type::element*)shared_mem;  // {blockDim.x, blockDim.y}

	// each thread loads one element from global to shared mem
	unsigned int ix  = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int tid = threadIdx.y;
	unsigned int iy  = blockIdx.y * blockDim.y + threadIdx.y;

	if(ix >= sizex)
		return;

	if(iy < sizey) {
		reduction_buffer[threadIdx.x + blockDim.x * tid] = kernel(ix, iy);
	} else {
		reduction_buffer[threadIdx.x + blockDim.x * tid] = (typename array_type::element)0.0;
	}

	__syncthreads();

	// do reduction in shared mem

	for(unsigned int s = blockDim.y / 2; s > 0; s >>= 1) {
		if(tid < s) {
			reduction_buffer[threadIdx.x + blockDim.x * tid] += reduction_buffer[threadIdx.x + blockDim.x * (tid + s)];
		}
		__syncthreads();
	}

	// write result for this block to global mem
	if(tid == 0)
		odata[blockIdx.y][ix] = reduction_buffer[threadIdx.x];
}

namespace gpu {
template<class type, size_t dim, class allocator = thrust::universal_allocator<type>>
using array = boost::multi::array<type, dim, allocator>;

struct reduce {
	explicit reduce(multi::ssize_t arg_size) : size(arg_size) {
	}
	multi::ssize_t size;
};

template<typename array_type>
struct array_access {
	array_type array;

	__host__ __device__ auto operator()(multi::ssize_t ii) const {
		return array[ii];
	}

	__host__ __device__ auto operator()(multi::ssize_t ix, multi::ssize_t iy) const {
		return array[ix][iy];
	}
};

template<class type, class kernel_type>
auto run(multi::ssize_t sizex, reduce const& redy, kernel_type kernel) /*-> gpu::array<decltype(kernel(0, 0)), 1>*/ {

	auto const sizey = redy.size;

	// using type = decltype(kernel(0, 0));

#ifndef ENABLE_GPU

	gpu::array<type, 1> accumulator(sizex, 0.0);

	for(multi::ssize_t iy = 0; iy < sizey; iy++) {
		for(multi::ssize_t ix = 0; ix < sizex; ix++) {
			accumulator[ix] += kernel(ix, iy);
		}
	}

	// thrust::transform(
	// 	boost::multi::extension_t({0, sizey}).begin(),
	// 	boost::multi::extension_t({0, sizey}).end(),
	// 	accumulator.begin(),
	// 	[](auto e) {return 0.0;}
	// );

	return accumulator;

#else

	gpu::array<type, 2> result;

	auto blocksize = 256;  // max_blocksize(reduce_kernel_vr<kernel_type, decltype(begin(result))>);

	unsigned bsizex = 4;  // this seems to be the optimal value
	if(sizex <= 2)
		bsizex = sizex;
	unsigned bsizey = blocksize / bsizex;

	assert(bsizey > 1);

	unsigned nblockx = (sizex + bsizex - 1) / bsizex;
	unsigned nblocky = (sizey + bsizey - 1) / bsizey;

	result.reextent({nblocky, sizex});

	gpu::array<type, 2> result2({nblocky, sizex});

	struct dim3 dg {
		nblockx, nblocky
	};
	struct dim3 db {
		bsizex, bsizey
	};

	auto shared_mem_size = blocksize * sizeof(type);

	assert(shared_mem_size <= 48 * 1024);

	reduce_kernel_vr<<<dg, db, shared_mem_size>>>(sizex, sizey, kernel, begin(result));
	// check_error(last_error());

	boost::multi::extents_t<2> xs = {sizex, sizey};
	assert( xs.size() == sizex );

	// thrust::transform_reduce(
	// 	xs.
	// );

	if(nblocky == 1) {
		cudaDeviceSynchronize();
		// gpu::sync();

		assert(result[0].size() == sizex);

		return gpu::array<type, 1>(result[0]);
	} else {
		return run<double>(sizex, reduce(nblocky), array_access<decltype(begin(result.transposed()))>{begin(result.transposed())});
	}

#endif
}

}  // namespace gpu

struct prod {
	/*__host__*/ __device__ auto operator()(multi::ssize_t ix, multi::ssize_t iy) const {
		return double(ix) * double(iy);
	}
};

auto main() -> int {
	#ifdef NDEBUG
	multi::ssize_t const maxsize = 39062;
	multi::ssize_t const nmax = 1000;
	#else
	multi::ssize_t const maxsize = 390;  // 390625;
	multi::ssize_t const nmax = 100;  // 10000;
	#endif
	auto pp = [] __host__ __device__(multi::index ix, multi::index iy) -> double { return double(ix) * double(iy); };

	std::chrono::microseconds mus{0};
	std::size_t FLOPs = 0;

	for(multi::ssize_t nx = 1; nx <= nmax; nx *= 10) {
		for(multi::ssize_t ny = 1; ny <= maxsize; ny *= 5) {

			multi::thrust::device_array<double, 2> M = [&]() {
				multi::thrust::universal_array<double, 2> ret({nx, ny});
				auto const [xs, ys] = ret.extents();
				for(auto ix : xs) {
					for(auto iy : ys) {
						ret[ix][iy] = pp(ix, iy);
					}
				}
				return ret;
			}();
			cudaDeviceSynchronize();
			DoNotOptimize(M);

			auto start = std::chrono::high_resolution_clock::now();

			multi::thrust::device_array<double, 1> sums(M.size());
			multi::thrust::reduce_by_index(M, sums);

			// auto sums = gpu::run<double>(nx, gpu::reduce(ny), gpu::array_access<decltype(M.begin())>{M.begin()});

			cudaDeviceSynchronize();
			DoNotOptimize(sums);

			mus += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
			FLOPs += nx*ny;
		}
	}
	std::cout << "time " << mus.count() << " µs FLOPS " << FLOPs << " flops/s " << FLOPs/(mus.count()/1e6)/1e9 << '\n';

	{
		auto nx = 13;
		auto ny = 19;

		multi::thrust::device_array<double, 1> sums(nx);
		multi::thrust::reduce_by_index(
			[] __host__ __device__(multi::index ix, multi::index iy) -> double { return double(ix) * double(iy); }
			^ multi::extents_t<2>(nx, ny),
			sums
		);
	}

	return boost::report_errors();
}

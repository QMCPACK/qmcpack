#ifdef COMPILATION // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
sudo cpupower frequency-set --governor performance; cmake -DCMAKE_CUDA_COMPILER=/usr/local/cuda-11.1/bin/nvcc -DCMAKE_BUILD_TYPE=Release .. && make $0.x && ctest ; exit
#endif // © Alfredo A. Correa 2021

#include <multi/array.hpp>
#include <multi/adaptors/blas/axpy.hpp>
#include <multi/adaptors/cuda/cublas/context.hpp>

#include <benchmark/benchmark.h>

#include <thrust/system/cuda/memory.h>

#include<execution>
#include<numeric>

namespace multi = boost::multi;

////////////////////////////////////////////////////////////////////////////////
template<class T, class X, class Y>
Y&& axpy_cpu_blas(T alpha, X const& x, Y&& y){
	return multi::blas::axpy(alpha, x, std::forward<Y>(y));
}

////////////////////////////////////////////////////////////////////////////////
template<class T, class X, class Y>
Y&& axpy_gpu_cublas(T alpha, X const& x, Y&& y){
	multi::cuda::cublas::context ctxt;
	return multi::blas::axpy(ctxt, alpha, x, std::forward<Y>(y));
}

////////////////////////////////////////////////////////////////////////////////
template<class T, class X, class Y>
Y&& axpy_cpu_transform(T alpha, X const& x, Y&& y){
	assert( size(x) == size(y) );
	std::transform( begin(x), end(x), begin(y), begin(y), [alpha](auto const& a, auto const& b){return a + alpha*b;} );
	return std::forward<Y>(y);
}

////////////////////////////////////////////////////////////////////////////////
template<class T> 
__global__ void axpy_gpu_legacy_kernel(int n, T alpha, T const* x, int x_stride, T* y, int y_stride){
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < n){*(y + i*y_stride) = alpha* *(x + i*x_stride) + *(y + i*y_stride);}
}

template<class T, class X, class Y>
Y&& axpy_gpu_legacy(T alpha, X const& x, Y&& y){
	assert( size(x) == size(y) );
	axpy_gpu_legacy_kernel<<<(size(x) + 255) / 256, 256>>>(size(x), alpha, raw_pointer_cast(base(x)), stride(x), raw_pointer_cast(base(y)), stride(y));
	return std::forward<Y>(y);
}

////////////////////////////////////////////////////////////////////////////////
template<class T, class ItX, class ItY> 
__global__ void axpy_gpu_multi_kernel(int n, T alpha, ItX x, ItY y){
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < n){y[i] = alpha*x[i] + y[i];}
}

template<class T, class X, class Y>
Y&& axpy_gpu_multi_kernel(T alpha, X const& x, Y&& y){
	assert( size(x) == size(y) );
	axpy_gpu_multi_kernel<<<(size(x) + 255) / 256, 256>>>(size(x), alpha, begin(x), begin(y));
	return std::forward<Y>(y);
}

////////////////////////////////////////////////////////////////////////////////

struct AXPY{
template<class T, class XV, class YV>
struct t{
	T alpha_;
	constexpr auto operator()(XV const& a, YV const& b) const{return alpha_*a + b;}
};

template<class T, class X, class Y, class XV = typename X::element_type, class YV = typename std::decay_t<Y>::element_type>
static Y&& gpu_multi_thrust(T alpha, X const& x, Y&& y){
	using thrust::transform;
	transform(begin(x), end(x), begin(y), begin(y), t<T, XV, YV>{alpha});//[alpha] (XV const& a, YV const& b){return alpha*a + b;} );
	return std::forward<Y>(y);
}
};

static void BM_axpy_cpu_transform(benchmark::State& state){

	using T = double;
	using alloc = std::allocator<double>;
	multi::array<T, 2, alloc> const X({1<<27, 2}, 1.);
	multi::array<T, 2, alloc>       Y(extensions(X), 1.);

	auto&& x = (~X)[0];
	auto&& y = (~Y)[0];

	for(auto _ : state){
		axpy_cpu_transform(2.0, x, y);
		benchmark::DoNotOptimize(base(y)); benchmark::ClobberMemory();
	}

	state.SetBytesProcessed(state.iterations()*size(x)*2.*sizeof(T));
	state.SetItemsProcessed(state.iterations()*size(x)             );
}
BENCHMARK(BM_axpy_cpu_transform);

static void BM_axpy_cpu_blas(benchmark::State& state){

	using T = double;
	using alloc = std::allocator<double>;
	multi::array<T, 2, alloc> const X({1<<27, 2}, 1.);
	multi::array<T, 2, alloc>       Y(extensions(X), 1.);

	auto&& x = (~X)[0];
	auto&& y = (~Y)[0];

	for(auto _ : state){
		axpy_cpu_blas(2.0, x, y);
		benchmark::DoNotOptimize(base(y)); benchmark::ClobberMemory();
	}

	state.SetBytesProcessed(state.iterations()*size(x)*2.*sizeof(T));
	state.SetItemsProcessed(state.iterations()*size(x)             );
}
BENCHMARK(BM_axpy_cpu_blas);

static void BM_axpy_gpu_cublas(benchmark::State& state){

	using T = double;
	using alloc = thrust::cuda::allocator<double>;
	multi::array<T, 2, alloc> const X({1<<27, 2}, 1.);
	multi::array<T, 2, alloc>       Y(extensions(X), 1.);

	auto&& x = (~X)[0];
	auto&& y = (~Y)[0];

	for(auto _ : state){
		axpy_gpu_cublas(2.0, x, y);
		benchmark::DoNotOptimize(base(y)); benchmark::ClobberMemory();
	}

	state.SetBytesProcessed(state.iterations()*size(x)*2.*sizeof(T));
	state.SetItemsProcessed(state.iterations()*size(x)             );
}
BENCHMARK(BM_axpy_gpu_cublas);

static void BM_axpy_gpu_legacy(benchmark::State& state){

	using T = double;
	using alloc = thrust::cuda::allocator<double>;
	multi::array<T, 2, alloc> const X({1<<27, 2}, 1.);
	multi::array<T, 2, alloc>       Y(extensions(X), 1.);

	auto&& x = (~X)[0];
	auto&& y = (~Y)[0];

	for(auto _ : state){
		axpy_gpu_legacy(2.0, x, y);
		cudaDeviceSynchronize(); benchmark::DoNotOptimize(y.base()); benchmark::ClobberMemory();
	}

	state.SetBytesProcessed(state.iterations()*size(x)*2.*sizeof(T));
	state.SetItemsProcessed(state.iterations()*size(x)             );
}
BENCHMARK(BM_axpy_gpu_legacy);

static void BM_axpy_gpu_multi_kernel(benchmark::State& state){

	using T = double;
	using alloc = thrust::cuda::allocator<double>;
	multi::array<T, 2, alloc> const X({1<<27, 2}, 1.);
	multi::array<T, 2, alloc>       Y(extensions(X), 1.);

	auto&& x = (~X)[0];
	auto&& y = (~Y)[0];

	for(auto _ : state){
		axpy_gpu_multi_kernel(2.0, x, y);
		cudaDeviceSynchronize(); benchmark::DoNotOptimize(base(y)); benchmark::ClobberMemory();
	}

	state.SetBytesProcessed(state.iterations()*size(x)*2.*sizeof(T));
	state.SetItemsProcessed(state.iterations()*size(x)             );
}
BENCHMARK(BM_axpy_gpu_multi_kernel);

static void BM_axpy_gpu_multi_thrust(benchmark::State& state){

	using T = double;
	using alloc = thrust::cuda::allocator<double>;
	multi::array<T, 2, alloc> const X({1<<27, 2}, 1.);
	multi::array<T, 2, alloc>       Y(extensions(X), 1.);

	auto&& x = (~X)[0];
	auto&& y = (~Y)[0];

	for(auto _ : state){
		AXPY::gpu_multi_thrust(2.0, x, y);
		cudaDeviceSynchronize(); benchmark::DoNotOptimize(base(y)); benchmark::ClobberMemory();
	}

	state.SetBytesProcessed(state.iterations()*size(x)*2.*sizeof(T));
	state.SetItemsProcessed(state.iterations()*size(x)             );
}
BENCHMARK(BM_axpy_gpu_multi_thrust);

BENCHMARK_MAIN();


#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
/usr/local/cuda-11.0/bin/nvcc -std=c++17 -ftemplate-backtrace-limit=0 $0 -o $0.$X `pkg-config --cflags --libs cudart-11.0 cuda-11.0` -lboost_timer&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2020

#include<benchmark/benchmark.h>

//#include<thrust/complex.h>
#include<thrust/device_allocator.h>
#include<thrust/device_vector.h>

#include "../../multi/array.hpp"
#include "../../multi/adaptors/thrust.hpp"

namespace multi = boost::multi;

#if not defined(NDEBUG)
#warning "Benchmark in debug mode?"
#endif

static void 
BM_cpu_vector_double_assignment
(benchmark::State& st){
	std::vector<double> const A(1<<28, 1.);
	std::vector<double>       B(A.size(), 2.);
	for(auto _ : st){
		B = A;
		benchmark::DoNotOptimize(B);
		benchmark::ClobberMemory();
	//	cudaDeviceSynchronize();
	}
	std::cout << A.size()*sizeof(A.front())/1e6 << "MB" << std::endl;
	st.SetBytesProcessed(st.iterations()*A.size()*sizeof(A.front()));
	st.SetItemsProcessed(st.iterations()*A.size());
}

static void 
BM_device_cudaMemcpy_double_assignment
(benchmark::State& st){
	thrust::device_vector<double> const A(1<<28, 1.);
	thrust::device_vector<double>       B(A.size(), 2.);
	for(auto _ : st){
		cudaMemcpy(raw_pointer_cast(B.data()), raw_pointer_cast(A.data()), A.size()*sizeof(A.front()), cudaMemcpyDeviceToDevice);
		cudaDeviceSynchronize();
	}
	st.SetBytesProcessed(st.iterations()*A.size()*sizeof(A.front()));
	st.SetItemsProcessed(st.iterations()*A.size());
}

static void 
BM_device_vector_double_assignment
(benchmark::State& st){
	thrust::device_vector<double> const A(1<<28, 1.);
	thrust::device_vector<double>       B(A.size(), 2.);
	for(auto _ : st){
		B = A;
	}
	st.SetBytesProcessed(st.iterations()*A.size()*sizeof(A.front()));
	st.SetItemsProcessed(st.iterations()*A.size());
}

static void BM_device_array_double_assignment(benchmark::State& st){
	using T = double;
	using alloc = thrust::device_allocator<T>; // std::allocator<T>;
	multi::array<T, 1, alloc> const A(1<<28, 1.);
	multi::array<T, 1, alloc>       B(extensions(A), 2.);
	for(auto _ : st){
		B() = A();
	}
	st.SetBytesProcessed(st.iterations()*A.num_elements()*sizeof(*A.base()));
	st.SetItemsProcessed(st.iterations()*A.num_elements());
}

static void BM_cpu_array_2D_double_assignment(benchmark::State& st){
	using T = double;
	using alloc = std::allocator<T>; //thrust::device_allocator<T>; // std::allocator<T>;
	multi::array<T, 2, alloc> const A({1<<14, 1<<14}, 1.);
	multi::array<T, 2, alloc>       B(extensions(A), 2.);
	for(auto _ : st){
		B() = A();
	}
	std::cout << A.num_elements()*sizeof(*A.base())/1e6 << "MB"<<std::endl;
	st.SetBytesProcessed(st.iterations()*A.num_elements()*sizeof(*A.base()));
	st.SetItemsProcessed(st.iterations()*A.num_elements());
}

static void BM_device_array_2D_double_assignment(benchmark::State& st){
	using T = double;
	using alloc = thrust::device_allocator<T>; // std::allocator<T>;
	multi::array<T, 2, alloc> const A({1<<14, 1<<14}, 1.);
	multi::array<T, 2, alloc>       B(extensions(A), 2.);
	for(auto _ : st){
		B() = A();
		cudaDeviceSynchronize();
	}
	if( B[10][10] == 2.) throw 0;
	st.SetBytesProcessed(st.iterations()*A.num_elements()*sizeof(*A.base()));
	st.SetItemsProcessed(st.iterations()*A.num_elements());
}

BENCHMARK(BM_cpu_vector_double_assignment);
BENCHMARK(BM_device_vector_double_assignment);
BENCHMARK(BM_device_cudaMemcpy_double_assignment);
BENCHMARK(BM_device_array_double_assignment);
BENCHMARK(BM_cpu_array_2D_double_assignment);
BENCHMARK(BM_device_array_2D_double_assignment);

BENCHMARK_MAIN();

//	using T = thrust::complex<double>;
//	using alloc = thrust::device_allocator<T>; // std::allocator<T>;
//	multi::array<T, 1, alloc> const A(10000, 1.);
//	multi::array<T, 1, alloc>       B(10000);
//	B = A;
//	assert( T{B[10]} == 1. );

//}


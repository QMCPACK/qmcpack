#ifdef COMPILATION// sudo cpupower frequency-set --governor performance
$CXX -DNDEBUG `#-DNOEXCEPT_ASSIGNMENT` $0 -o $0x `pkg-config --libs benchmark`&&$0x&&rm $0x;exit
#endif

#include <benchmark/benchmark.h>

#include "../array.hpp"

#include<complex>
#include<numeric>

namespace multi = boost::multi;

using complex = std::complex<double>; 
MAYBE_UNUSED constexpr complex I{0, 1};
using T = complex;

constexpr std::size_t N = 1 << 24;

static void VectorAssignment(benchmark::State& state){
	std::vector<T> a(N); std::iota(begin(a), end(a), 1.11);
	std::vector<T> b(N); std::iota(begin(b), end(b), 0.);
	for(auto _ : state){
		b = a;
		benchmark::DoNotOptimize(b.data());
		benchmark::ClobberMemory();
	}
}
BENCHMARK(VectorAssignment)->Unit(benchmark::kMillisecond);

static void MultiAssignment(benchmark::State& state){
	multi::array<T, 1> a(N); std::iota(begin(a), end(a), 1.11);
	multi::array<T, 1> b(N); std::iota(begin(b), end(b), 0.);
	for(auto _ : state){
		b = a;
		benchmark::DoNotOptimize(b.data());
		benchmark::ClobberMemory();
	}
}
BENCHMARK(MultiAssignment)->Unit(benchmark::kMillisecond);

static void RawAssignment(benchmark::State& state){
	multi::array<T, 1> a(N); std::iota(begin(a), end(a), 1.11);
	multi::array<T, 1> b(N); std::iota(begin(b), end(b), 0.);
	for(auto _ : state){
		std::copy_n(a.data_elements(), a.num_elements(), b.data_elements());
		benchmark::DoNotOptimize(b.data());
		benchmark::ClobberMemory();
	}
}
BENCHMARK(RawAssignment)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

//2020-06-12 12:01:19
//Running ./assignment.cppx
//Run on (12 X 4600 MHz CPU s)
//CPU Caches:
//  L1 Data 32K (x6)
//  L1 Instruction 32K (x6)
//  L2 Unified 256K (x6)
//  L3 Unified 12288K (x1)
//Load Average: 12.13, 9.36, 7.01
//-----------------------------------------------------------
//Benchmark                 Time             CPU   Iterations
//-----------------------------------------------------------
//VectorAssignment       50.0 ms         26.7 ms           25
//MultiAssignment        53.7 ms         28.5 ms           26
//RawAssignment          52.0 ms         27.8 ms           23


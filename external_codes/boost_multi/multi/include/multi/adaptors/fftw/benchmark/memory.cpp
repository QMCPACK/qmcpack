#ifdef COMPILATION// sudo cpupower frequency-set --governor performance && sudo apt install libbenchmark-dev
${CXX:-c++} -Ofast -DNDEBUG -march=native `#-DNOEXCEPT_ASSIGNMENT` -I../../../../../include/ $0 -o $0x `pkg-config --libs benchmark fftw3`&&$0x&&rm $0x;exit
#endif

#include <benchmark/benchmark.h>

#include <multi/array.hpp>
#include <multi/adaptors/fftw/memory.hpp>
#include <multi/adaptors/fftw.hpp>

#include<complex>

namespace multi = boost::multi;

using complex = std::complex<double>;

template <class Alloc>
static void Allocation(benchmark::State& state){

    multi::array<complex, 2, Alloc> in({state.range(0), state.range(0)*2}, 1.2);
	multi::array<complex, 2, Alloc> out(extensions(in), 3.1);

    std::vector<double> v(state.range(0)*3.14);
    benchmark::DoNotOptimize(v);

    benchmark::DoNotOptimize(in);
    benchmark::DoNotOptimize(out);

    benchmark::ClobberMemory();

    multi::fftw::plan p(std::array<bool, 2>{true, true}, in, out, multi::fftw::forward, multi::fftw::estimate);
	for(auto _ : state){
		benchmark::DoNotOptimize(in);
		benchmark::DoNotOptimize(out);
		// benchmark::ClobberMemory();

        p();
	}

    benchmark::DoNotOptimize(in);
    benchmark::DoNotOptimize(out);
    benchmark::ClobberMemory();
}

BENCHMARK(Allocation<std        ::allocator<complex>>)->DenseRange(100, 500, 28);  // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
BENCHMARK(Allocation<multi::fftw::allocator<complex>>)->DenseRange(100, 500, 28);  // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

template <class Alloc>
static void Allocation1D(benchmark::State& state){

    multi::array<complex, 1, Alloc> in({state.range(0)}, 1.2);
	multi::array<complex, 1, Alloc> out(extensions(in), 3.1);

    std::vector<double> v(state.range(0)*3.14);
    benchmark::DoNotOptimize(v);

    benchmark::DoNotOptimize(in);
    benchmark::DoNotOptimize(out);

    benchmark::ClobberMemory();

    multi::fftw::plan p(std::array<bool, 1>{true}, in, out, multi::fftw::forward, multi::fftw::estimate);
	for(auto _ : state){
		benchmark::DoNotOptimize(in);
		benchmark::DoNotOptimize(out);
		// benchmark::ClobberMemory();

        p();
	}

    benchmark::DoNotOptimize(in);
    benchmark::DoNotOptimize(out);
    benchmark::ClobberMemory();
}

BENCHMARK(Allocation1D<std        ::allocator<complex>>)->RangeMultiplier(2)->Range(128, 128*1024);
BENCHMARK(Allocation1D<multi::fftw::allocator<complex>>)->RangeMultiplier(2)->Range(128, 128*1024);

BENCHMARK_MAIN();

/*
$ sh ./memory.cpp
2023-09-08T10:11:42-07:00
Running ./memory.cppx
Run on (12 X 4000.06 MHz CPU s)
CPU Caches:
  L1 Data 32 KiB (x6)
  L1 Instruction 32 KiB (x6)
  L2 Unified 256 KiB (x6)
  L3 Unified 12288 KiB (x1)
Load Average: 4.28, 3.36, 2.76
***WARNING*** Library was built as DEBUG. Timings may be affected.
-----------------------------------------------------------------------------------------------
Benchmark                                                     Time             CPU   Iterations
-----------------------------------------------------------------------------------------------
Allocation<std ::allocator<complex>>/100                 106167 ns       105362 ns         6801
Allocation<std ::allocator<complex>>/128                 186345 ns       186233 ns         3701
Allocation<std ::allocator<complex>>/156                 445409 ns       445387 ns         1609
Allocation<std ::allocator<complex>>/184                1524822 ns      1524758 ns          423
Allocation<std ::allocator<complex>>/212                6334445 ns      6334253 ns          116
Allocation<std ::allocator<complex>>/240                1639634 ns      1639581 ns          454
Allocation<std ::allocator<complex>>/268                7807228 ns      7806998 ns           97
Allocation<std ::allocator<complex>>/296               13261983 ns     13261620 ns           45
Allocation<std ::allocator<complex>>/324                3183764 ns      3183686 ns          200
Allocation<std ::allocator<complex>>/352                4003407 ns      4003311 ns          174
Allocation<std ::allocator<complex>>/380                7525147 ns      7524610 ns          106
Allocation<std ::allocator<complex>>/408               10952047 ns     10951038 ns           51
Allocation<std ::allocator<complex>>/436               42780039 ns     42778256 ns           16
Allocation<std ::allocator<complex>>/464               19106268 ns     19104952 ns           36
Allocation<std ::allocator<complex>>/492               50939459 ns     50933888 ns           13
Allocation<multi::fftw::allocator<complex>>/100          161494 ns       161399 ns         4491
Allocation<multi::fftw::allocator<complex>>/128          242198 ns       242146 ns         2874
Allocation<multi::fftw::allocator<complex>>/156          661061 ns       660975 ns         1024
Allocation<multi::fftw::allocator<complex>>/184         2139026 ns      2138641 ns          322
Allocation<multi::fftw::allocator<complex>>/212         9557704 ns      9556300 ns           77
Allocation<multi::fftw::allocator<complex>>/240         2313975 ns      2313203 ns          285
Allocation<multi::fftw::allocator<complex>>/268        11401940 ns     11399821 ns           64
Allocation<multi::fftw::allocator<complex>>/296        19879090 ns     19874211 ns           35
Allocation<multi::fftw::allocator<complex>>/324         3994153 ns      3990909 ns          176
Allocation<multi::fftw::allocator<complex>>/352         5923172 ns      5923003 ns          120
Allocation<multi::fftw::allocator<complex>>/380         9400210 ns      9398525 ns           73
Allocation<multi::fftw::allocator<complex>>/408         9982523 ns      9981809 ns           67
Allocation<multi::fftw::allocator<complex>>/436        42858298 ns     42854934 ns           17
Allocation<multi::fftw::allocator<complex>>/464        20109321 ns     20104686 ns           32
Allocation<multi::fftw::allocator<complex>>/492        54813030 ns     54097910 ns           13
Allocation1D<std ::allocator<complex>>/128                  307 ns          307 ns      2482984
Allocation1D<std ::allocator<complex>>/256                  597 ns          597 ns      1167097
Allocation1D<std ::allocator<complex>>/512                 1325 ns         1325 ns       446234
Allocation1D<std ::allocator<complex>>/1024                2833 ns         2833 ns       226483
Allocation1D<std ::allocator<complex>>/2048                7819 ns         7817 ns        88305
Allocation1D<std ::allocator<complex>>/4096               23896 ns        23894 ns        27796
Allocation1D<std ::allocator<complex>>/8192               64528 ns        64520 ns        10695
Allocation1D<std ::allocator<complex>>/16384             122521 ns       122508 ns         5928
Allocation1D<std ::allocator<complex>>/32768             294487 ns       294375 ns         2577
Allocation1D<std ::allocator<complex>>/65536             635923 ns       635821 ns         1090
Allocation1D<std ::allocator<complex>>/131072           1399999 ns      1399604 ns          486
Allocation1D<multi::fftw::allocator<complex>>/128           285 ns          285 ns      2578003
Allocation1D<multi::fftw::allocator<complex>>/256           604 ns          604 ns      1140959
Allocation1D<multi::fftw::allocator<complex>>/512          1206 ns         1206 ns       473758
Allocation1D<multi::fftw::allocator<complex>>/1024         2761 ns         2759 ns       246249
Allocation1D<multi::fftw::allocator<complex>>/2048         7240 ns         7235 ns        89839
Allocation1D<multi::fftw::allocator<complex>>/4096        20633 ns        20611 ns        31992
Allocation1D<multi::fftw::allocator<complex>>/8192        56934 ns        56813 ns        12723
Allocation1D<multi::fftw::allocator<complex>>/16384      121458 ns       121206 ns         5358
Allocation1D<multi::fftw::allocator<complex>>/32768      281074 ns       280545 ns         2503
Allocation1D<multi::fftw::allocator<complex>>/65536      606937 ns       606822 ns         1098
Allocation1D<multi::fftw::allocator<complex>>/131072    1582752 ns      1582109 ns          478
*/

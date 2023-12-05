#ifdef COMPILATION_INSTRUCTIONS
clang++-9 -std=c++17 -O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS $0 -o $0x -lbenchmark && $0x $@ && rm $0x; exit
#endif

#include <boost/multi_array.hpp>

#include <benchmark/benchmark.h>

#include "../include/multi/array.hpp"

const long X_SIZE = 4000;
const long Y_SIZE = 4000;

typedef boost::multi_array<double, 2> ImageArrayType;

static void MeasureNative2(benchmark::State& state) {

    // Create the native array
    double *nativeMatrix = new double [X_SIZE * Y_SIZE]{};

	for (auto _ : state) {
        for (long y = 0; y != Y_SIZE; ++y)
        {
            for (long x = 0; x != X_SIZE; ++x)
            {
                nativeMatrix[x + (y * X_SIZE)] *= 2.345;
            }
        }
		benchmark::DoNotOptimize(nativeMatrix);
	}
	delete[] nativeMatrix;
}
BENCHMARK(MeasureNative2);

static void MeasureNative(benchmark::State& state) {

    // Create the native array
    double *nativeMatrix = new double [X_SIZE * Y_SIZE]{};

	for (auto _ : state) {
        for (long y = 0; y != Y_SIZE; ++y)
        {
            for (long x = 0; x != X_SIZE; ++x)
            {
                nativeMatrix[x + (y * X_SIZE)] *= 2.345;
            }
        }
		benchmark::DoNotOptimize(nativeMatrix);
	}
	delete[] nativeMatrix;
}
BENCHMARK(MeasureNative);

static void MeasureBoostMultiArray(benchmark::State& state) {

    ImageArrayType boostMatrix(boost::extents[X_SIZE][Y_SIZE]);

	for (auto _ : state) {
        for (long y = 0; y != Y_SIZE; ++y)
        {
            for (long x = 0; x != X_SIZE; ++x)
            {
                boostMatrix[x][y] *= 2.345;
            }
        }
		benchmark::DoNotOptimize(boostMatrix);
	}
}
BENCHMARK(MeasureBoostMultiArray);


static void MeasureBoostMultiArrayInverted(benchmark::State& state) {

    ImageArrayType boostMatrix(boost::extents[X_SIZE][Y_SIZE]);

	for (auto _ : state) {
        for (long x = 0; x != X_SIZE; ++x)
        {
            for (long y = 0; y != Y_SIZE; ++y)
            {
                boostMatrix[x][y] *= 2.345;
            }
        }
		benchmark::DoNotOptimize(boostMatrix);
	}
}
BENCHMARK(MeasureBoostMultiArrayInverted);

static void MeasureBoostMultiArrayRaw(benchmark::State& state) {

    ImageArrayType boostMatrix(boost::extents[X_SIZE][Y_SIZE]);

	for (auto _ : state) {
        for (unsigned long n = 0; n != boostMatrix.num_elements(); ++n)
        {
            boostMatrix.data()[n] *= 2.345;
        }
		benchmark::DoNotOptimize(boostMatrix);
	}
}
BENCHMARK(MeasureBoostMultiArrayRaw);

namespace multi = boost::multi;

static void MeasureMulti(benchmark::State& state) {

    multi::array<double, 2> M({X_SIZE, Y_SIZE}, 0.);

	for (auto _ : state) {
        for (long y = 0; y != Y_SIZE; ++y)
        {
            for (long x = 0; x != X_SIZE; ++x)
            {
                M[x][y] *= 2.345;
            }
        }
		benchmark::DoNotOptimize(M);
	}
}
BENCHMARK(MeasureMulti);

static void MeasureMultiInverted(benchmark::State& state) {

    multi::array<double, 2> M({X_SIZE, Y_SIZE}, 0.);

	for (auto _ : state) {
        for (long x = 0; x != X_SIZE; ++x)
        {
            for (long y = 0; y != Y_SIZE; ++y)
            {
                M[x][y] *= 2.345;
            }
        }
		benchmark::DoNotOptimize(M);
	}
}
BENCHMARK(MeasureMultiInverted);

static void MeasureMultiRaw(benchmark::State& state) {

    multi::array<double, 2> M({X_SIZE, Y_SIZE}, 0.);

	for (auto _ : state) {
        for (unsigned long n = 0; n != M.num_elements(); ++n)
        {
            M.data_elements()[n] *= 2.345;
        }
		benchmark::DoNotOptimize(M);
	}
}
BENCHMARK(MeasureBoostMultiArrayRaw);

BENCHMARK_MAIN();


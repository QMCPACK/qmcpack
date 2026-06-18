// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// © Alfredo A. Correa 2019-2021

#include <benchmark/benchmark.h>

#include "../array.hpp"
#include "../algorithms/gemm.hpp"
#include "../adaptors/blas/gemm.hpp"

namespace multi = boost::multi;

template<class MatrixA, class MatrixB, class MatrixC>
auto naive_product(MatrixA const& A, MatrixB const& B, double beta, MatrixC&& C) -> MatrixC&& {
	assert(   C .size() ==   A. size() );
	assert( (~C).size() == (~B).size() );
	for(auto i : extension(C)) {
		for(auto j : extension(~C)) {
			C[i][j] = std::inner_product(begin(A[i]), end(A[i]), begin((~B)[j]), beta*C[i][j]);
		}
	}
	return std::forward<MatrixC>(C);
}

auto const N = 1024;

static void Bnaive_product(benchmark::State& _) {
	auto const A = []{
		multi::array<double, 2> A({       N,        N}    );
		std::iota(begin(A.elements()), end(A.elements()), 0.1);
		return A;
	}();

	auto const B = [&]{
		multi::array<double, 2> B({(~A).size(),        N}    );
		std::iota(begin(B.elements()), end(B.elements()), 0.2);
		return B;
	}();

	multi::array<double, 2> C({  A.size() , (~B).size()}, 0.);

    benchmark::ClobberMemory();
	while(_.KeepRunning()) {
		naive_product(A, B, 0., C);
		benchmark::DoNotOptimize(C.data_elements());
	    benchmark::ClobberMemory();
	}
}

static void multi_detail_naive_gemm(benchmark::State& state) {
	auto const A = []{
		multi::array<double, 2> A({       N,        N}    );
		std::iota(begin(A.elements()), end(A.elements()), 0.1);
		return A;
	}();

	auto const B = [&]{
		multi::array<double, 2> B({(~A).size(),        N}    );
		std::iota(begin(B.elements()), end(B.elements()), 0.2);
		return B;
	}();

	multi::array<double, 2> C({  A.size() , (~B).size()}, 0.);

    benchmark::ClobberMemory();
	for(auto _ : state) {
		multi::detail::naive_gemm(1., A, B, 0., C);
		benchmark::DoNotOptimize(C.data_elements());
	    benchmark::ClobberMemory();
	}
}

static void multi_gemm(benchmark::State& state) {
	auto const A = []{
		multi::array<double, 2> A({       N,        N}    );
		std::iota(begin(A.elements()), end(A.elements()), 0.1);
		return A;
	}();

	auto const B = [&]{
		multi::array<double, 2> B({(~A).size(),        N}    );
		std::iota(begin(B.elements()), end(B.elements()), 0.2);
		return B;
	}();

	multi::array<double, 2> C({  A.size() , (~B).size()}, 0.);

    benchmark::ClobberMemory();
	for(auto _ : state) {
		multi::gemm(1., A, B, 0., C);
		benchmark::DoNotOptimize(C);
	    benchmark::ClobberMemory();
	}
}

static void multi_blas_gemm(benchmark::State& state) {
	auto const A = [] {
		multi::array<double, 2> A({       N,        N}    );
		std::iota(begin(A.elements()), end(A.elements()), 0.1);
		return A;
	}();

	auto const B = [&] {
		multi::array<double, 2> B({(~A).size(),        N}    );
		std::iota(begin(B.elements()), end(B.elements()), 0.2);
		return B;
	}();

	multi::array<double, 2> C({A.size() , (~B).size()}, 0.);

    benchmark::ClobberMemory();
	for(auto _ : state) {
		multi::blas::gemm(1., A, B, 0., C);
		benchmark::DoNotOptimize(C);
	    benchmark::ClobberMemory();
	}
}

BENCHMARK(Bnaive_product);

BENCHMARK(multi_detail_naive_gemm);
BENCHMARK(multi_gemm);

BENCHMARK(multi_blas_gemm);

BENCHMARK_MAIN();

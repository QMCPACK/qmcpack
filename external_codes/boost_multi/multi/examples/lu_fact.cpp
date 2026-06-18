#ifdef COMPILATION  // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 - std = c++ 17 - o $0x - lboost_timer `pkg - config-- libs tbb` && $0x && rm $0x;
exit
#endif
// Copyright 2018-2024 Alfredo A. Correa

#include <boost/multi/array.hpp>

#include <boost/timer/timer.hpp>

#include <algorithm>  // transform
#include <iostream>
#include <numeric>  // iota
#include <tuple>

namespace multi = boost::multi;

template<class Matrix>
Matrix&& lu_fact(Matrix&& A) {
	using multi::size;
	auto m = size(A);  // n = size(A[0]);//std::get<1>(sizes(A));
	using multi::sizes;
	using std::begin;
	using std::end;
	for(auto k = 0 * m; k != std::min(m - 1, get<1>(sizes(A))); ++k) {
		auto const& Ak  = A[k];
		auto const& Akk = Ak[k];
		std::for_each(
			begin(A) + k + 1, end(A), [&](auto&& Ai) {
				std::transform(
					begin(Ai) + k + 1, end(Ai), begin(Ak) + k + 1, begin(Ai) + k + 1,
					[z = (Ai[k] /= Akk)](auto a, auto b) { return a - z * b; }
				);
			}
		);
	}
	return std::forward<Matrix>(A);
}

template<class Matrix>
Matrix&& lu_fact2(Matrix&& A) {
	using multi::size;
	auto const [m, n] = A.sizes();

	for(decltype(m) k = 0; k != m - 1; ++k) {
		for(auto i = k + 1; i != m; ++i){
			auto const z = A[i][k]/A[k][k];
			A[i][k] = z;
			std::transform(begin(A[i]) + k + 1, begin(A[i]) + std::max(n, k + 1), A[k].begin() + k + 1, begin(A[i]) + k + 1, [&](auto a, auto b) { return a - z * b; });
		}
	}
	return std::forward<Matrix>(A);
}

template<class Matrix>
Matrix&& lu_fact3(Matrix&& A) {
	using multi::size;
	auto const [m, n] = A.sizes();
	for(auto k = 0 * m; k != m - 1; ++k) {
		auto&& Ak = A[k];
		std::for_each(
			begin(A) + k + 1, end(A), [&](auto& Ai) {
				auto const z = Ai[k] / Ak[k];
				Ai[k]        = z;
				assert((k + 1) <= n);
				for(auto j = k + 1; j < n; ++j)
					Ai[j] -= z * Ak[j];
			}
		);
	}
	return std::forward<Matrix>(A);
}

using std::cout;
int main() {
	{
		multi::array<double, 2> A = {
			{-3.0, 2.0, -4.0},
			{ 0.0, 1.0,  2.0},
			{ 2.0, 4.0,  5.0},
		};
		multi::array<double, 1> y = {12.0, 5.0, 2.0};
		double AA[3][3];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy types
		using std::copy;
		copy(begin(A), end(A), begin(*multi::array_ptr(&AA)));

		lu_fact(A);
		lu_fact(AA);
		assert(std::equal(begin(A), end(A), begin(*multi::array_ptr(&AA)), end(*multi::array_ptr(&AA))));
	}
	{
		multi::array<double, 2> A({6000, 7000});
		std::iota(A.elements().begin(), A.elements().end(), 0.1);
		std::transform(A.elements().begin(), A.elements().begin() + A.num_elements(), A.elements().begin(), [](auto x) { return x /= 2.0e6; });
		{
			boost::timer::auto_cpu_timer t;
			lu_fact(A({3000, 6000}, {0, 4000}));
			cout << A[456][123] << std::endl;
		}
	}
}

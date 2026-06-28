#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX -DNDEBUG $0 -o $0x&&$0x&&rm $0x;exit
#endif
// Copyright 2019-2023 Alfredo A. Correa

#include <multi/array.hpp>

#include<algorithm>
#include<chrono>
#include<iostream>
#include<numeric> // iota
#include<string>
#include<utility> // move
#include<vector>

namespace multi = boost::multi;
using std::cout;

template<class Matrix, class Vector, class idx = typename std::decay_t<Vector>::difference_type>
auto gj_solve(Matrix&& A, Vector&& y) -> decltype(y[0] /= A[0][0], y) {
	idx Asize = size(A);
	for(idx r = 0; r != Asize; ++r) {
		auto&&      Ar  = A[r];
		auto const& Arr = Ar[r];
		for(idx c = r + 1; c != Asize; ++c)
			Ar[c] /= Arr;
		auto const& yr = (y[r] /= Arr);
		for(idx r2 = r + 1; r2 != Asize; ++r2) {
			auto&&      Ar2  = A[r2];
			auto const& Ar2r = A[r2][r];
			for(idx c = r + 1; c != Asize; ++c)
				Ar2[c] -= Ar2r * Ar[c];
			y[r2] -= Ar2r * yr;
		}
	}
	for(idx r = Asize - 1; r > 0; --r) {
		auto const& yr = y[r];
		for(idx r2 = r - 1; r2 >= 0; --r2)
			y[r2] -= A[r2][r] * yr;
	}
	return y;
}

template<class Matrix, class Vector, class idx = typename std::decay_t<Vector>::difference_type>
auto gj_solve2(Matrix&& A, Vector&& y) -> decltype(y[0] /= A[0][0], y) {
	idx Asize = size(A);
	for(idx r = 0; r != Asize; ++r) {
		auto&&      Ar  = A[r];
		auto const& Arr = Ar[r];
		//  std::transform(Ar.begin() + r + 1, Ar.end(), Ar.begin() + r + 1, [&](auto const& a){return a/Arr;});
		for(idx c = r + 1; c != Asize; ++c)
			Ar[c] /= Arr;
		auto const& yr = (y[r] /= Arr);
		for(idx r2 = r + 1; r2 != Asize; ++r2) {
			auto&&      Ar2  = A[r2];
			auto const& Ar2r = A[r2][r];
			std::transform(std::move(Ar2).begin() + r + 1, std::move(Ar2).end(), std::move(Ar).begin() + r + 1, std::move(Ar2).begin() + r + 1, [&](auto&& a, auto&& b) { return a - Ar2r * b; });
			y[r2] -= Ar2r * yr;
		}
	}
	for(idx r = Asize - 1; r > 0; --r) {
		auto const& yr = y[r];
		for(idx r2 = r - 1; r2 >= 0; --r2)
			y[r2] -= A[r2][r] * yr;
	}
	return y;
}

namespace {
// minimal RAII wall-clock timer (replaces boost::timer::auto_cpu_timer)
class auto_timer {
	std::string                           label_;
	std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();

 public:
	explicit auto_timer(std::string label = {}) : label_{std::move(label)} {}
	auto_timer(auto_timer const&)                    = delete;
	auto_timer(auto_timer&&)                         = delete;
	auto operator=(auto_timer const&) -> auto_timer& = delete;
	auto operator=(auto_timer&&) -> auto_timer&      = delete;
	~auto_timer() {
		std::cerr << label_ << std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count() << " s (wall)\n";
	}
};
}  // namespace

int main(){
	{
		multi::array<double, 2> A = {
			{-3.0, 2.0, -4.0},
			{ 0.0, 1.0,  2.0},
			{ 2.0, 4.0,  5.0},
		};
		multi::array<double, 1> y = {12.0, 5.0, 2.0};  //(M); assert(y.size() == M); iota(y.begin(), y.end(), 3.1);
		gj_solve(A, y);
		cout << y[0] << " " << y[1] << " " << y[2] << std::endl;
	}
	{
		multi::array<double, 2> A({6000, 7000});
		std::iota(A.data_elements(), A.data_elements() + A.num_elements(), 0.1);
		std::transform(A.data_elements(), A.data_elements() + A.num_elements(), A.data_elements(), [](auto x) { return x /= 2.e6; });
		std::vector<double> y(3000);
		std::iota(y.begin(), y.end(), 0.2);
		{
			auto_timer t;
			gj_solve(A({1000, 4000}, {0, 3000}), y);
		}
		cout << y[45] << std::endl;
	}
	{
		multi::array<double, 2> A({6000, 7000});
		std::iota(A.data_elements(), A.data_elements() + A.num_elements(), 0.1);
		std::transform(A.data_elements(), A.data_elements() + A.num_elements(), A.data_elements(), [](auto x) { return x /= 2.e6; });
		std::vector<double> y(3000);
		std::iota(y.begin(), y.end(), 0.2);
		{
			auto_timer t;
			gj_solve2(A({1000, 4000}, {0, 3000}), y);
		}
		cout << y[45] << std::endl;
	}
}


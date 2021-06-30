#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX -DNDEBUG $0 -o $0x -lboost_timer&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#include "../array.hpp"

#include<iostream>
#include<vector>
#include<numeric> // iota
#include<algorithm>

namespace multi = boost::multi;
using std::cout;

template<class Matrix, class Vector, class idx = typename std::decay_t<Vector>::difference_type>
auto gj_solve(Matrix&& A, Vector&& y)->decltype(y[0]/=A[0][0], y){
	idx Asize = size(A);
	for(idx r = 0; r != Asize; ++r){
		auto&& Ar = A[r]; auto const& Arr = Ar[r];
		for(idx c = r + 1; c != Asize; ++c) Ar[c] /= Arr;
		auto const& yr = (y[r] /= Arr);
		for(idx r2 = r + 1; r2 != Asize; ++r2){
			auto&& Ar2 = A[r2]; auto const& Ar2r = A[r2][r];
			for(idx c = r + 1; c != Asize; ++c) Ar2[c] -= Ar2r*Ar[c];
			y[r2] -= Ar2r*yr;
		}
	}
	for(idx r = Asize - 1; r > 0; --r){
		auto const& yr = y[r];
		for(idx r2 = r-1; r2 >=0; --r2) y[r2] -= A[r2][r]*yr;
	}
	return y;
}

template<class Matrix, class Vector, class idx = typename std::decay_t<Vector>::difference_type>
auto gj_solve2(Matrix&& A, Vector&& y)->decltype(y[0]/=A[0][0], y){
	idx Asize = size(A);
	for(idx r = 0; r != Asize; ++r){
		auto&& Ar = A[r]; auto const& Arr = Ar[r];
	//	std::transform(Ar.begin() + r + 1, Ar.end(), Ar.begin() + r + 1, [&](auto const& a){return a/Arr;});
		for(idx c = r + 1; c != Asize; ++c) Ar[c] /= Arr;
		auto const& yr = (y[r] /= Arr);
		for(idx r2 = r + 1; r2 != Asize; ++r2){
			auto&& Ar2 = A[r2]; auto const& Ar2r = A[r2][r];
			std::transform(std::move(Ar2).begin() + r + 1, std::move(Ar2).end(), std::move(Ar).begin() + r + 1, std::move(Ar2).begin() + r + 1, [&](auto&& a, auto&& b){return a - Ar2r*b;});
			y[r2] -= Ar2r*yr;
		}
	}
	for(idx r = Asize - 1; r > 0; --r){
		auto const& yr = y[r];
		for(idx r2 = r-1; r2 >=0; --r2) y[r2] -= A[r2][r]*yr;
	}
	return y;
}

#include <boost/timer/timer.hpp>

int main(){
	{
		multi::array<double, 2> A = {{-3., 2., -4.},{0., 1., 2.},{2., 4., 5.}};
		multi::array<double, 1> y = {12.,5.,2.}; //(M); assert(y.size() == M); iota(y.begin(), y.end(), 3.1);
		gj_solve(A, y);
		cout << y[0] <<" "<< y[1] <<" "<< y[2] << std::endl;
	}
	{
		multi::array<double, 2> A({6000, 7000}); std::iota(A.data(), A.data() + A.num_elements(), 0.1);
		std::transform(A.data(), A.data() + A.num_elements(), A.data(), [](auto x){return x/=2.e6;});
		std::vector<double> y(3000); std::iota(y.begin(), y.end(), 0.2);
		{
			boost::timer::auto_cpu_timer t;
			gj_solve(A({1000, 4000}, {0, 3000}), y);
		}
		cout << y[45] << std::endl;
	}
	{
		multi::array<double, 2> A({6000, 7000}); std::iota(A.data(), A.data() + A.num_elements(), 0.1);
		std::transform(A.data(), A.data() + A.num_elements(), A.data(), [](auto x){return x/=2.e6;});
		std::vector<double> y(3000); std::iota(y.begin(), y.end(), 0.2);
		{
			boost::timer::auto_cpu_timer t;
			gj_solve2(A({1000, 4000}, {0, 3000}), y);
		}
		cout << y[45] << std::endl;
	}
}


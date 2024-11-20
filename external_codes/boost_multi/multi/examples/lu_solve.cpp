#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX -DNDEBUG $0 -o $0x -lboost_timer&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2021

#include "../array.hpp"

#include<cmath> // for std::abs
#include<numeric> // for std::iota

struct lup{ // LU method for decomposition and solution

// translated from  https://en.wikipedia.org/wiki/LU_decomposition#C_code_example
template<class Matrix, class Permutation>
static auto decompose(Matrix&& A, Permutation&& P, double tol = std::numeric_limits<double>::epsilon()){
	std::iota(begin(P), end(P), typename std::decay_t<Permutation>::value_type{0});
	auto const N = std::min(size(A), size(~A));
	assert( P.size() >= N );

	auto&& ret = A({0, N}, {0, N});
	for(auto i : extension(ret)){
		if(lup::permute_max_diagonal(A, P, i) < tol) return A({0, i}, {0, i});

		for(auto&& row : A({i + 1, N})){
			auto&& urow = row({i + 1, N});
			std::transform(
				cbegin(urow), cend(urow), cbegin(A[i]({i + 1, N})), begin(urow), 
				[f = row[i] /= A[i][i]](auto const& a, auto const& b){return a - f*b;}
			);
		}
	}
	return std::move(ret);
}

template<class Matrix, class Permutation, class VectorSol>
static auto solve(Matrix const& LU, Permutation const& P, VectorSol&& x) -> VectorSol&&{
	return upper_solve(LU, lower_solve(LU, permute(P, x)));
}

private:

template<class Matrix, class Permutation, class Index>
static auto permute_max_diagonal(Matrix&& LU, Permutation&& P, Index i){
	auto mi = std::max_element(begin(LU) + i, end(LU), [i](auto const& a, auto const& b){return std::abs(a[i]) < std::abs(b[i]);}) - begin(LU);
	     swap(LU[i], LU[mi]); 
	std::swap(P [i], P [mi]);
	return std::abs(LU[i][i]);
}

template<class Permutation, class Vector>
static auto permute(Permutation const& p, Vector&& data) -> Vector&&{
	assert(size(p) <= size(data));
	using index = typename Permutation::size_type;
	for(index i = 0; i != size(p); ++i){
		index k = p[i];
		for( ; k > i; k = p[k]){}
		index pk = p[k];
		if(k >=i and pk != i){
			auto const t = data[i];
			for( ; pk != i; k = pk, pk = p[k]){
				data[k] = data[pk];
			};
			data[k] = t;
		}
	}
	return std::forward<Vector>(data);
}

template<class LUMatrix, class Vector>
static auto lower_solve(LUMatrix const& LU, Vector&& x) -> Vector&&{
	assert(size(LU) <= size(x));
	auto const N = size(LU);
	for(typename LUMatrix::size_type i = 0; i != N; ++i){
		auto const& Lrowi = LU[i]({0, i});
		x[i] -= std::inner_product(begin(Lrowi), end(Lrowi), cbegin(x), 0.);
	}
	return std::forward<Vector>(x);
}

template<class LUMatrix, class Vector>
static auto upper_solve(LUMatrix const& LU, Vector&& x) -> Vector&&{
	assert(size(LU) <= size(x));
	auto const N = size(LU);
	for(typename LUMatrix::size_type i = N - 1; i >= 0; --i){
		auto const& Urowi = LU[i]({i + 1, N});
		(x[i] -= std::inner_product(begin(Urowi), end(Urowi), cbegin(x) + i + 1, 0.)) /= LU[i][i];
	}
	return std::forward<Vector>(x);
}

};

namespace multi = boost::multi;

int main(){
	multi::array<double, 2> const Aconst = {
		{ 6.80, -6.05, -0.45, 8.32, -9.67},
		{-2.11, -3.30,  2.58, 2.71, -5.14},
		{ 5.66,  5.36, -2.70, 4.35, -7.26},
		{ 5.97, -4.44,  0.27,-7.17,  6.08},
		{ 8.23,  1.08,  9.04, 2.14, -6.87}
	};
	auto A = Aconst;
	multi::array<int, 1> P({5}, 0.);
	lup::decompose(A, P);

	multi::array<double, 1> x = {4.02, 6.19, -8.22, -7.57, -3.03};

	lup::solve(A, P, x);
	
	assert( std::abs(x[4] - 0.565756) < 1e-4);
}


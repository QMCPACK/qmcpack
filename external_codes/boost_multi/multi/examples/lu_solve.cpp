#ifdef COMPILATION  // clang-format off
set -x;${CXX:-c++} -std=c++20 -O3 -DNDEBUG -I../include $0 -o $0x &&time $0x&&rm $0x;exit
#endif  // clang-format on

// Copyright 2021-2024 Alfredo A. Correa

#include <boost/multi/array.hpp>

#include <cmath>  // for std::abs
#include <iostream>
#include <numeric>  // for std::iota
#include <random>
#include <source_location>

	namespace {
	struct lup {  // LU method for decomposition and solution

		// translated from  https://en.wikipedia.org/wiki/LU_decomposition#C_code_example
		template<class Matrix, class Permutation>
		static auto decompose(Matrix&& A, Permutation&& P, double tol = std::numeric_limits<double>::epsilon()) {
			auto const N = std::min(size(A), size(~A));
			assert(P.size() >= N);

			if(A.size() == 0)
				return A({0, 1}, {0, 1});

			auto&& ret = A({0, N}, {0, N});

			// if(lup::permute_max(A, P) < tol) return A({0, 1}, {0, 1});
			// for(auto&& row : A({1, N})) {
			//  auto&& urow = row({1, N});
			//  std::transform(
			//      cbegin(urow), cend(urow), cbegin(A[0]({1, N})), begin(urow),
			//      [f = row[i] /= A[0][0]](auto const& a, auto const& b) { return a - f * b; }
			//  );
			// }

			// decompose(A({1, N}, {1, N}), P({1, N}));

			for(auto i : extension(ret)) {
				auto&& A_rest = A({i, N});
				auto&& P_rest = P({i, N});
				auto&& A_rest_rest = A_rest({0, N - i}, {i, N});
				if(lup::permute_max_diagonal(A_rest({0, N - i}, {0, N}), P_rest({0, N - i}), i) < tol)
					return A({0, i}, {0, i});

				for(auto&& row : A_rest({1, N - i})) {
					auto&& urow = row({i + 1, N});
					std::transform(
						cbegin(urow), cend(urow), cbegin(A_rest_rest[0]({1, N - i})), begin(urow),
						[f = row[i] /= A_rest_rest[0][0]](auto const& a, auto const& b) { return a - f * b; }
					);
				}
			}
			return std::move(ret);
		}

		// translated from  https://en.wikipedia.org/wiki/LU_decomposition#C_code_example
		template<
			class Permutation, class Matrix>
		static auto decompose(Matrix&& A, double tol = std::numeric_limits<double>::epsilon()) {
			Permutation P(A.size());
			std::iota(begin(P), end(P), typename std::decay_t<Permutation>::value_type{0});

			decompose(A, P, tol);
			return std::pair<Matrix&&, Permutation>(std::forward<Matrix>(A), std::move(P));
		}

		template<class Matrix, class Permutation, class VectorSol>
		static auto solve(Matrix const& LU, Permutation const& P, VectorSol&& x) -> VectorSol&& {
			return upper_solve(LU, lower_solve(LU, permute(P, x)));
		}

	 private:
		template<class Matrix, class Permutation, class Index>
		static auto permute_max_diagonal(Matrix&& LU, Permutation&& P, Index i) {
			auto mi = std::max_element(begin(LU), end(LU), [i](auto const& a, auto const& b) { return std::abs(a[i]) < std::abs(b[i]); }) - begin(LU);
			using std::swap;
			swap(LU[0], LU[mi]);
			swap(P[0], P[mi]);
			return std::abs(LU[0][i]);
		}

		template<class Matrix, class Permutation, class Index>
		static auto permute_max(Matrix&& LU, Permutation&& P) {
			auto mi = std::max_element(begin(LU), end(LU), [](auto const& a, auto const& b) { return std::abs(a[0]) < std::abs(b[0]); }) - begin(LU);
			if(mi != 0) {
				using std::swap;
				swap(LU[0], LU[mi]);
				swap(P[0], P[mi]);
			}
			return std::abs(LU[0][0]);
		}

		template<class Permutation, class Vector>
		static auto permute(Permutation const& p, Vector&& data) -> Vector&& {
			assert(size(p) <= size(data));
			using index = typename Permutation::size_type;
			for(index i = 0; i != size(p); ++i) {
				index k = p[i];
				for(; k > i; k = p[k]) {
				}
				index pk = p[k];
				if(k >= i and pk != i) {
					auto const t = data[i];
					for(; pk != i; k = pk, pk = p[k]) {
						data[k] = data[pk];
					};
					data[k] = t;
				}
			}
			return std::forward<Vector>(data);
		}

		template<class LUMatrix, class Vector>
		static auto lower_solve(LUMatrix const& LU, Vector&& x) -> Vector&& {
			assert(size(LU) <= size(x));
			auto const N = size(LU);
			for(typename LUMatrix::size_type i = 0; i != N; ++i) {
				auto const& Lrowi = LU[i]({0, i});
				x[i] -= std::inner_product(begin(Lrowi), end(Lrowi), cbegin(x), 0.);
			}
			return std::forward<Vector>(x);
		}

		template<class LUMatrix, class Vector>
		static auto upper_solve(LUMatrix const& LU, Vector&& x) -> Vector&& {
			assert(size(LU) <= size(x));
			auto const N = size(LU);
			for(typename LUMatrix::size_type i = N - 1; i >= 0; --i) {
				auto const& Urowi = LU[i]({i + 1, N});
				(x[i] -= std::inner_product(begin(Urowi), end(Urowi), cbegin(x) + i + 1, 0.)) /= LU[i][i];
			}
			return std::forward<Vector>(x);
		}
	};
}

namespace multi = boost::multi;

int main() try {
	{
		multi::array<double, 2> const Aconst = {
			{ 6.80, -6.05, -0.45,  8.32, -9.67},
			{-2.11, -3.30,  2.58,  2.71, -5.14},
			{ 5.66,  5.36, -2.70,  4.35, -7.26},
			{ 5.97, -4.44,  0.27, -7.17,  6.08},
			{ 8.23,  1.08,  9.04,  2.14, -6.87}
		};
		auto A = Aconst;
		// multi::array<int, 1> P(multi::array<int, 1>::extensions_type{5});
		auto P = std::get<1>(lup::decompose<multi::array<int, 1>>(A));

		multi::array<double, 1> x = {4.02, 6.19, -8.22, -7.57, -3.03};

		lup::solve(A, P, x);
	
		std::cout << std::abs(x[4]) << std::endl;
		(std::abs(x[4] - 0.565756) < 1e-4) ?: throw std::source_location::current();
	}
	{
		multi::array<double, 2> A({4000, 4000});

		std::random_device rd;         // Non-deterministic random number generator
		std::mt19937       gen(rd());  // Mersenne Twister engine seeded with rd()

		std::uniform_real_distribution<double> dist(-1.0, 1.0);
		for(auto&& e : A.elements()) {
			e = dist(gen);
		}

		multi::array<int, 1> P(multi::array<int, 1>::extensions_type{A.size()});
		std::iota(P.begin(), P.end(), int{});
		lup::decompose(A, P);

		multi::array<double, 1> x(A.size());
		for(auto&& e : x.elements()) {
			e = dist(gen);
		}

		lup::solve(A, P, x);
	}
} catch(std::source_location const& loc) {
	std::cerr << loc.file_name() << ':' << loc.line() << '\n';
	throw;
}

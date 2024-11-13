// Copyright 2019-2024 Alfredo A. Correa

#include <boost/multi/adaptors/lapack/filling.hpp>  // for filling, filling...
#include <boost/multi/adaptors/lapack/potrf.hpp>    // for potrf

// IWYU pragma: no_include "boost/multi/adaptors/blas/complex_traits.hpp"  // for blas
#include <boost/multi/adaptors/blas/gemm.hpp>        // for gemm
#include <boost/multi/adaptors/blas/herk.hpp>        // for herk
// IWYU pragma: no_include "boost/multi/adaptors/blas/numeric.hpp"     // for underlying
#include <boost/multi/adaptors/blas/operations.hpp>  // for H, (anonymous)

#include <boost/multi/array.hpp>                     // for array, subarray

#include <algorithm>  // for for_each, generate  // IWYU pragma: keep
#include <complex>    // for operator*, complex
#include <cstdlib>    // for abs
#include <iostream>   // for operator<<, ostream
#include <limits>     // for numeric_limits
#include <random>     // for uniform_real_dis...
#include <string>     // for allocator, opera...
// IWYU pragma: no_include <tuple>
// IWYU pragma: no_include <type_traits>  // for add_const<>::type
#include <utility>  // for forward

namespace multi = boost::multi;
// namespace lapack = multi::lapack;
namespace blas = multi::blas;

using complex = std::complex<double>;

auto operator<<(std::ostream& os, std::complex<double> const& cx) -> std::ostream& {
	return os << real(cx) << " + I*" << imag(cx);
}

template<class M> auto print(M const& arr) -> decltype(auto) { return print(arr, ""); }
template<class M> auto print(M const& arr, std::string const& msg) -> decltype(auto) {
	using multi::size;
	using std::cout;
	cout << msg << "\n"
		 << '{';
	for(int i = 0; i != size(arr); ++i) {
		cout << '{';
		for(auto j : arr[i].extension()) {  // NOLINT(altera-unroll-loops)
			cout << arr[i][j];
			if(j + 1 != size(arr[i])) {
				cout << ", ";
			}
		}
		cout << '}' << '\n';
		if(i + 1 != size(arr)) {
			cout << ", ";
		}
	}
	return cout << '}' << '\n';
}

template<class M>
auto print(M const& arr, char const* msg) -> decltype(auto) {
	return print(arr, std::string{msg});  // NOLINT(fuchsia-default-arguments-calls)
}

template<class M>
auto randomize(M&& arr) -> M&& {
	std::random_device dev;
	std::mt19937       eng{dev()};

	auto gen = [&]() {
		auto unif = std::uniform_real_distribution<>{-1.0, 1.0};
		return std::complex<double>(unif(eng), unif(eng));
	};

	std::for_each(begin(arr), end(arr), [&](auto&& row) { std::generate(begin(row), end(row), gen); });
	return std::forward<M>(arr);
}

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/
#define BOOST_TEST_CLOSE(X, Y, ToL) BOOST_TEST(std::abs((X) - (Y)) < (ToL))

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	/*
	BOOST_AUTO_TEST_CASE(orthogonalization_over_rows, *boost::unit_test::tolerance(0.00001)){
		auto A = randomize(multi::array<complex, 2>({3, 10}));
		lapack::onrm(A);

		using blas::herk;
		using blas::hermitized;
		using blas::filling;
		auto id = herk(filling::upper, A);
		BOOST_TEST( real(id[1][1]) == 1.0 ); BOOST_TEST( imag(id[1][1]) == 0.0 );
		BOOST_TEST( real(id[1][2]) == 0.0 ); BOOST_TEST( imag(id[1][2]) == 0.0 );
	}
	*/

	// BOOST_AUTO_TEST_CASE(orthogonalization_over_rows_cuda, *boost::unit_test::tolerance(0.00001)) {
	//  auto Acpu = randomize(multi::array<complex, 2>({3, 10}));

	//  multi::cuda::array<complex, 2> A = Acpu;

	//  using namespace blas;
	//  using namespace lapack;

	//  trsm(filling::lower, hermitized(potrf(filling::upper, herk(filling::upper, A))), A);

	//  Acpu    = A;
	//  auto id = herk(filling::upper, Acpu);
	//  BOOST_TEST( real(id[1][1]) == 1.0 );
	//  BOOST_TEST( imag(id[1][1]) == 0.0 );
	//  BOOST_TEST( real(id[1][2]) == 0.0 );
	//  BOOST_TEST( imag(id[1][2]) == 0.0 );
	// }

	/*
	BOOST_AUTO_TEST_CASE(orthogonalization_over_columns, *boost::unit_test::tolerance(0.00001)){

		auto A = randomize( multi::array<complex, 2>({10, 3}) );
		using blas::hermitized;
		lapack::onrm(hermitized(A));

		using blas::filling;
		auto id = herk(filling::upper, hermitized(A));
		BOOST_TEST( real(id[1][1]) == 1. ); BOOST_TEST( imag(id[1][1]) == 0. );
		BOOST_TEST( real(id[1][2]) == 0. ); BOOST_TEST( imag(id[1][2]) == 0. );
	}*/

	BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_define_both_sides) {  // }, *boost::unit_test::tolerance(0.0000001)) {
		auto const   I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

		multi::array<complex, 2> const A_gold = {
			{3.23 + 0.00 * I,  1.51 - 1.92 * I,  1.90 + 0.84 * I,  0.42 + 2.50 * I},
			{1.51 + 1.92 * I,  3.58 + 0.00 * I, -0.23 + 1.11 * I, -1.18 + 1.37 * I},
			{1.90 - 0.84 * I, -0.23 - 1.11 * I,  4.09 + 0.00 * I,  2.33 - 0.14 * I},
			{0.42 - 2.50 * I, -1.18 - 1.37 * I,  2.33 + 0.14 * I,  4.29 + 0.00 * I},
		};

		auto A = A_gold;  // NOLINT(readability-identifier-length) lapack conventional name

		auto const As = multi::lapack::potrf(multi::lapack::filling::upper, A).size();
		BOOST_TEST( As == A.size() );

		auto AA = A;

		for(auto i = 0; i != 4; ++i) {
			for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
				AA[i][j] = 0.0;
			}
		}

		auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

		for(auto i = 0; i != 4; ++i) {
			for(auto j = 0; j != 4; ++j) {  // NOLINT(altera-unroll-loops)
				BOOST_TEST_CLOSE(real(A_gold[i][j]), real(C[i][j]), 0.0000001);
				BOOST_TEST_CLOSE(imag(A_gold[i][j]), imag(C[i][j]), 0.0000001);
			}
		}
	}

	BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_define_upper) {
		double const nan = std::numeric_limits<double>::quiet_NaN();
		auto const   I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

		multi::array<complex, 2> const A_gold = {
			{3.23 + 0.00 * I, 1.51 - 1.92 * I,  1.90 + 0.84 * I,  0.42 + 2.50 * I},
			{  nan + nan * I, 3.58 + 0.00 * I, -0.23 + 1.11 * I, -1.18 + 1.37 * I},
			{  nan - nan * I,   nan - nan * I,  4.09 + 0.00 * I,  2.33 - 0.14 * I},
			{  nan - nan * I,   nan - nan * I,    nan + nan * I,  4.29 + 0.00 * I},
		};

		auto A = A_gold;  // NOLINT(readability-identifier-length) lapack conventional name

		auto const As = multi::lapack::potrf(multi::lapack::filling::upper, A).size();

		BOOST_TEST( As == A.size() );

		auto AA = A;

		for(auto i = 0; i != 4; ++i) {
			for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
				AA[i][j] = 0.0;
			}
		}

		auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

		print(A_gold, "A gold");  // NOLINT(fuchsia-default-arguments-calls)
		print(C, "recover");      // NOLINT(fuchsia-default-arguments-calls)

		for(auto i = 0; i != 4; ++i) {
			// only compare upper part of the reference array (the other half is garbage)
			for(auto j = i; j != 4; ++j) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)
				BOOST_TEST_CLOSE(real(A_gold[i][j]), real(C[i][j]), 0.0000001);
				// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay,readability-simplify-boolean-expr) bug in clang-tidy 14
				BOOST_TEST_CLOSE(imag(A_gold[i][j]), imag(C[i][j]), 0.0000001);
			}
		}
	}

	BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_trivial_imperfect) {
		double const nan = std::numeric_limits<double>::quiet_NaN();
		auto const   I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

		multi::array<complex, 2> const A_gold = {
			{3.23 + 0.00 * I, 1.51 - 1.92 * I,      1.90 + 0.84 * I,     0.42 + 2.50 * I},
			{  nan + nan * I, 3.58 + 0.00 * I,     -0.23 + 1.11 * I,    -1.18 + 1.37 * I},
			{  nan - nan * I,   nan - nan * I, -10000.00 + 0.00 * I,     0.00 - 0.00 * I},
			{  nan - nan * I,   nan - nan * I,        nan + nan * I, -1000.00 + 0.00 * I},
		};

		auto A = A_gold;  // NOLINT(readability-identifier-length) lapack conventional name

		auto const& Adec = multi::lapack::potrf(multi::lapack::filling::upper, A);

		print(A, "A");
		print(Adec, "A dec");

		auto AA = +Adec;

		// NOLINTNEXTLINE(altera-id-dependent-backward-branch)
		for(auto i = 0; i != AA.size(); ++i) {
			for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
				AA[i][j] = 0.0;
			}
		}

		auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

		print(A_gold, "A gold");  // NOLINT(fuchsia-default-arguments-calls)
		print(C, "recover");      // NOLINT(fuchsia-default-arguments-calls)

		using std::get;

		for(auto i = 0; i != AA.size(); ++i) {  // NOLINT(altera-id-dependent-backward-branch)
			for(auto j = i; j != get<1>(C.sizes()); ++j) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)
				BOOST_TEST_CLOSE(real(A_gold[i][j]), real(C[i][j]), 0.0000001);
				// NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-array-to-pointer-decay,readability-simplify-boolean-exp) bug in clang-tidy 14
				BOOST_TEST_CLOSE(imag(A_gold[i][j]), imag(C[i][j]), 0.0000001);
			}
		}
	}

	BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_nontrivial_imperfect) {
		double const nan = std::numeric_limits<double>::quiet_NaN();
		auto const   I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

		multi::array<complex, 2> const A_gold = {
			{1.00 + 0.00 * I, 0.00 - 0.00 * I,  0.00 + 0.00 * I,  0.00 + 0.00 * I},
			{  nan + nan * I, 1.00 + 0.00 * I,  0.00 + 0.00 * I,  0.00 + 0.00 * I},
			{  nan - nan * I,   nan - nan * I, -1.00 + 0.00 * I,  0.00 - 0.00 * I},
			{  nan - nan * I,   nan - nan * I,    nan + nan * I, -1.00 + 0.00 * I},
		};

		auto A = A_gold;  // NOLINT(readability-identifier-length) lapack conventional name

		auto const& Adec = multi::lapack::potrf(multi::lapack::filling::upper, A);

		print(A, "A");
		print(Adec, "A dec");

		auto AA = +Adec;  // NOLINT(altera-id-dependent-backward-branch) bug in clang-tidy 14

		// NOLINTNEXTLINE(altera-id-dependent-backward-branch)
		for(auto i = 0; i != AA.size(); ++i) {
			for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
				AA[i][j] = 0.0;
			}
		}

		auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

		print(A_gold, "A gold");  // NOLINT(fuchsia-default-arguments-calls)
		print(C, "recover");      // NOLINT(fuchsia-default-arguments-calls)

		using std::get;  // workaround use of function template name with no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		// NOLINTNEXTLINE(altera-id-dependent-backward-branch)
		for(auto i = 0; i != AA.size(); ++i) {
			// NOLINTNEXTLINE(altera-unroll-loops,altera-id-dependent-backward-branch)
			for(auto j = i; j != get<1>(C.sizes()); ++j) {  // only compare upper part of the reference array (the other half is garbage)
				BOOST_TEST_CLOSE(real(A_gold[i][j]), real(C[i][j]), 0.0000001);
				BOOST_TEST_CLOSE(imag(A_gold[i][j]), imag(C[i][j]), 0.0000001);
			}
		}
	}

	BOOST_AUTO_TEST_CASE(lapack_potrf) {  // , *boost::unit_test::tolerance(0.00001)) {
		double const nan = std::numeric_limits<double>::quiet_NaN();
		auto const   I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

		{
			// NOLINTNEXTLINE(readability-identifier-length)
			multi::array<complex, 2> A = {
				{167.413 + 0.0 * I, 126.804 - 0.00143505 * I, 125.114 - 0.1485590 * I},
				{    nan + nan * I,        167.381 + 0.0 * I, 126.746 + 0.0327519 * I},
				{    nan + nan * I,            nan + nan * I,       167.231 + 0.0 * I},
			};

			print(A, "original A");
			using boost::multi::lapack::filling;
			using boost::multi::lapack::potrf;

			auto const As = potrf(filling::upper, A).size();  // A is hermitic in upper triangular (implicit below)
			BOOST_TEST( As == A.size() );

			BOOST_TEST_CLOSE(real(A[1][2]), 3.78646, 0.00001);
			BOOST_TEST_CLOSE(imag(A[1][2]), 0.0170734, 0.00001);
			//  BOOST_TEST( A[2][1] != A[2][1] );
			print(A, "decomposition");

			multi::array<complex, 2> C(A.extensions(), complex{0.0, 0.0});  // NOLINT(readability-identifier-length) conventional lapack name

			multi::array<complex, 2> AA = A;

			auto const [is, js] = AA.extensions();
			for(auto i : is) {
				for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch)
					AA[i][j] = std::conj(A[j][i]);
				}
			}

			blas::gemm(complex{1.0, 0.0}, blas::H(AA), AA, complex{0.0, 0.0}, C);

			print(C, "recovery");
		}
		// {
		//  multi::cuda::managed::array<complex, 2> A = {
		//      {167.413, 126.804 - 0.00143505 * I, 125.114 - 0.1485590 * I},
		//      {    NAN,                  167.381, 126.746 + 0.0327519 * I},
		//      {    NAN,                      NAN,                 167.231},
		//  };
		//  using lapack::filling;
		//  using lapack::potrf;
		//  potrf(filling::upper, A);  // A is hermitic in upper triangular (implicit below)
		//  BOOST_TEST( real(A[1][2]) == 3.78646 );
		//  BOOST_TEST( imag(A[1][2]) == 0.0170734 );
		//  //  BOOST_TEST( A[2][1] != A[2][1] );
		// }
		// {
		//  multi::cuda::array<complex, 2> A = {
		//      {167.413, 126.804 - 0.00143505 * I, 125.114 - 0.1485590 * I},
		//      {    NAN,                  167.381, 126.746 + 0.0327519 * I},
		//      {    NAN,                      NAN,                 167.231},
		//  };
		//  using lapack::filling;
		//  using lapack::potrf;
		//  potrf(filling::upper, A);  // A is hermitic in upper triangular (implicit below)
		//  multi::array<complex, 2> A_copy = A;
		//  print(A_copy);
		// }
	}

	return boost::report_errors();
}

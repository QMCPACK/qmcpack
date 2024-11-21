// Copyright 2019-2024 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuSolver potrf"
#include <boost/test/unit_test.hpp>

#include "../../blas/gemm.hpp"
#include "../../blas/herk.hpp"

#include "../../lapack/potrf.hpp"

#include <iostream>
#include <random>

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
template<class M> auto print(M const& arr, char const* msg) -> decltype(auto) { return print(arr, std::string{msg}); }  // NOLINT(fuchsia-default-arguments-calls)

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

BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_define_both_sides, *boost::unit_test::tolerance(0.0000001)) {
	double const nan = std::numeric_limits<double>::quiet_NaN();
	auto const   I   = complex{0.0, 1.0};  // NOLINT(readability-identifier-length)

	multi::array<complex, 2> const A_gold = {
		{3.23 + 0.00 * I,  1.51 - 1.92 * I,  1.90 + 0.84 * I,  0.42 + 2.50 * I},
		{1.51 + 1.92 * I,  3.58 + 0.00 * I, -0.23 + 1.11 * I, -1.18 + 1.37 * I},
		{1.90 - 0.84 * I, -0.23 - 1.11 * I,  4.09 + 0.00 * I,  2.33 - 0.14 * I},
		{0.42 - 2.50 * I, -1.18 - 1.37 * I,  2.33 + 0.14 * I,  4.29 + 0.00 * I},
	};

	auto A = A_gold;  // NOLINT(readability-identifier-length) lapack conventional name

	auto const As = multi::lapack::potrf(multi::lapack::filling::upper, A).size();
	BOOST_REQUIRE( As == A.size() );

	auto AA = A;

	for(auto i = 0; i != 4; ++i) {
		for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
			AA[i][j] = 0.0;
		}
	}

	auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

	for(auto i = 0; i != 4; ++i) {
		for(auto j = 0; j != 4; ++j) {
			BOOST_TEST( real(A_gold[i][j]) == real(C[i][j]) );
			BOOST_TEST( imag(A_gold[i][j]) == imag(C[i][j]) );
		}
	}
}

BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_define_upper, *boost::unit_test::tolerance(0.0000001)) {
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

	BOOST_REQUIRE( As == A.size() );

	auto AA = A;

	for(auto i = 0; i != 4; ++i) {
		for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
			AA[i][j] = 0.0;
		}
	}

	auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

	print(A_gold, "A gold");  // NOLINT(fuchsia-default-arguments-calls)
	print(C, "recover");  // NOLINT(fuchsia-default-arguments-calls)

	for(auto i = 0; i != 4; ++i) {
		for(auto j = i; j != 4; ++j) {  // NOLINT(altera-id-dependent-backward-branch)  // only compare upper part of the reference array (the other half is garbage)
			BOOST_TEST( real(A_gold[i][j]) == real(C[i][j]) );
			BOOST_TEST( imag(A_gold[i][j]) == imag(C[i][j]) );
		}
	}
}

BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_trivial_imperfect, *boost::unit_test::tolerance(0.0000001)) {  // NOLINT(fuchsia-default-arguments-calls)
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

	for(auto i = 0; i != AA.size(); ++i) {  // NOLINT(altera-id-dependent-backward-branch)
		for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
			AA[i][j] = 0.0;
		}
	}

	auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

	print(A_gold, "A gold");  // NOLINT(fuchsia-default-arguments-calls)
	print(C, "recover");  // NOLINT(fuchsia-default-arguments-calls)

	for(auto i = 0; i != AA.size(); ++i) {  // NOLINT(altera-id-dependent-backward-branch)
		for(auto j = i; j != std::get<1>(C.sizes()); ++j) {  // only compare upper part of the reference array (the other half is garbage)  // NOLINT(altera-id-dependent-backward-branch)
			BOOST_TEST( real(A_gold[i][j]) == real(C[i][j]) );
			BOOST_TEST( imag(A_gold[i][j]) == imag(C[i][j]) );
		}
	}
}

BOOST_AUTO_TEST_CASE(numericalalgorithmsgroup_nontrivial_imperfect, *boost::unit_test::tolerance(0.0000001)) {  // NOLINT(fuchsia-default-arguments-calls)
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

	auto AA = +Adec;

	for(auto i = 0; i != AA.size(); ++i) {  // NOLINT(altera-id-dependent-backward-branch)
		for(auto j = 0; j != i; ++j) {  // NOLINT(altera-unroll-loops)
			AA[i][j] = 0.0;
		}
	}

	auto const C = +blas::herk(1.0, blas::H(AA));  // +blas::gemm(1.0, blas::H(AA), AA);  // NOLINT(readability-identifier-length) conventional lapack name

	print(A_gold, "A gold");  // NOLINT(fuchsia-default-arguments-calls)
	print(C, "recover");  // NOLINT(fuchsia-default-arguments-calls)

	for(auto i = 0; i != AA.size(); ++i) {  // NOLINT(altera-id-dependent-backward-branch)
		for(auto j = i; j != std::get<1>(C.sizes()); ++j) {  // only compare upper part of the reference array (the other half is garbage)  // NOLINT(altera-id-dependent-backward-branch)
			BOOST_TEST( real(A_gold[i][j]) == real(C[i][j]) );
			BOOST_TEST( imag(A_gold[i][j]) == imag(C[i][j]) );
		}
	}
}

BOOST_AUTO_TEST_CASE(lapack_potrf, *boost::unit_test::tolerance(0.00001)) {
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
		BOOST_REQUIRE( As == A.size() );

		BOOST_TEST( real(A[1][2]) == 3.78646 );
		BOOST_TEST( imag(A[1][2]) == 0.0170734 );
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

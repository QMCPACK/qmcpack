#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&nvcc -x cu --expt-relaxed-constexpr`#$CXX` $0 -o $0x -Wno-deprecated-declarations -lcudart -lcublas -lcusolver `pkg-config --libs blas lapack` -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework -DBOOST_LOG_DYN_LINK -lboost_log -lpthread -lboost_system &&$0x&&rm $0x $0.cpp; exit
#endif
// © Alfredo A. Correa 2019-2020
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuSolver potrf"
#include<boost/test/unit_test.hpp>

#include "../../../adaptors/cuda.hpp"
#include             "../../lapack/potrf.hpp"
#include             "../../blas/herk.hpp"
#include             "../../blas/trsm.hpp"

#include "../../../adaptors/lapack/cuda.hpp"
#include "../../../adaptors/blas/cuda.hpp"

#include<iostream>
#include<random>

namespace multi = boost::multi;
namespace lapack = multi::lapack;
namespace blas = multi::blas;

using complex = std::complex<double>;

std::ostream& operator<<(std::ostream& os, std::complex<double> const& c){
	return os<< real(c) <<" + I*"<< imag(c);
}

template<class M> decltype(auto) print(M const& C){
	using std::cout;
	using multi::size;
	cout<<'{';
	for(int i = 0; i != size(C); ++i){
		cout<<'{';
		for(int j = 0; j != size(C[i]); ++j){
			cout<< C[i][j];
			if(j + 1 != size(C[i])) cout<<", ";
		}
		cout<<'}'<<std::endl;
		if(i + 1 != size(C)) cout<<", ";
	}
	return cout<<'}'<<std::endl;
}

template<class M>
M&& randomize(M&& A){
	std::mt19937 eng{123};
	auto gen = [&](){return std::complex<double>{std::uniform_real_distribution<>{-1, 1}(eng), std::uniform_real_distribution<>{-1, 1}(eng)};};
	std::for_each(begin(A), end(A), [&](auto&& r){std::generate(begin(r), end(r), gen);});
	return std::forward<M>(A);
}

/*
BOOST_AUTO_TEST_CASE(orthogonalization_over_rows, *boost::unit_test::tolerance(0.00001)){
	auto A = randomize(multi::array<complex, 2>({3, 10}));
	lapack::onrm(A);

	using blas::herk;
	using blas::hermitized;
	using blas::filling;
	auto id = herk(filling::upper, A);
	BOOST_TEST( real(id[1][1]) == 1. ); BOOST_TEST( imag(id[1][1]) == 0. );
	BOOST_TEST( real(id[1][2]) == 0. ); BOOST_TEST( imag(id[1][2]) == 0. );
}
*/

BOOST_AUTO_TEST_CASE(orthogonalization_over_rows_cuda, *boost::unit_test::tolerance(0.00001)){
	auto Acpu = randomize(multi::array<complex, 2>({3, 10}));

	multi::cuda::array<complex, 2> A = Acpu;

	using namespace blas;
	using namespace lapack;

	trsm(filling::lower, hermitized(potrf(filling::upper, herk(filling::upper, A))), A);

	Acpu = A;
	auto id = herk(filling::upper, Acpu);
	BOOST_TEST( real(id[1][1]) == 1. ); BOOST_TEST( imag(id[1][1]) == 0. );
	BOOST_TEST( real(id[1][2]) == 0. ); BOOST_TEST( imag(id[1][2]) == 0. );
}

/*
BOOST_AUTO_TEST_CASE(orthogonalization_over_columns, *boost::unit_test::tolerance(0.00001)){

	auto A = randomize(	multi::array<complex, 2>({10, 3}) );
	using blas::hermitized;
	lapack::onrm(hermitized(A));

	using blas::filling;
	auto id = herk(filling::upper, hermitized(A));
	BOOST_TEST( real(id[1][1]) == 1. ); BOOST_TEST( imag(id[1][1]) == 0. );
	BOOST_TEST( real(id[1][2]) == 0. ); BOOST_TEST( imag(id[1][2]) == 0. );
}*/

BOOST_AUTO_TEST_CASE(lapack_potrf, *boost::unit_test::tolerance(0.00001) ){

	complex const I{0, 1};
{
	multi::array<complex, 2> A = {
		{167.413, 126.804 - 0.00143505*I, 125.114 - 0.1485590*I},
		{NAN    , 167.381               , 126.746 + 0.0327519*I},
		{NAN    , NAN                   , 167.231              }
	};
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::upper, A); // A is hermitic in upper triangular (implicit below)
	BOOST_TEST( real(A[1][2]) == 3.78646 );
	BOOST_TEST( imag(A[1][2]) == 0.0170734 );
//	BOOST_TEST( A[2][1] != A[2][1] );
	print(A);
}
{
	multi::cuda::managed::array<complex, 2> A = {
		{167.413, 126.804 - 0.00143505*I, 125.114 - 0.1485590*I},
		{NAN    , 167.381               , 126.746 + 0.0327519*I},
		{NAN    , NAN                   , 167.231              }
	};
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::upper, A); // A is hermitic in upper triangular (implicit below)
	BOOST_TEST( real(A[1][2]) == 3.78646 );
	BOOST_TEST( imag(A[1][2]) == 0.0170734 );
//	BOOST_TEST( A[2][1] != A[2][1] );
}
{
	multi::cuda::array<complex, 2> A = {
		{167.413, 126.804 - 0.00143505*I, 125.114 - 0.1485590*I},
		{NAN    , 167.381               , 126.746 + 0.0327519*I},
		{NAN    , NAN                   , 167.231              }
	};
	using lapack::filling;
	using lapack::potrf;
	potrf(filling::upper, A); // A is hermitic in upper triangular (implicit below)
	multi::array<complex, 2> A_copy = A;
	print(A_copy);
}

}


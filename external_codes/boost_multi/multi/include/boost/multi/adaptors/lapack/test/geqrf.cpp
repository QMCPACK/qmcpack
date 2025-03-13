#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&nvcc -x cu --expt-relaxed-constexpr`#$CXX` $0 -o $0x -Wno-deprecated-declarations -lcudart -lcublas -lcusolver `pkg-config --libs blas lapack` -DBOOST_TEST_DYN_LINK -lboost_unit_test_framework -DBOOST_LOG_DYN_LINK -lboost_log -lpthread -lboost_system &&$0x&&rm $0x $0.cpp; exit
#endif
// © Alfredo A. Correa 2019-2020
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi getrf"
#include<boost/test/unit_test.hpp>

#include<multi/array.hpp>
#include<multi/adaptors/lapack/getrf.hpp>
#include<multi/adaptors/blas/gemm.hpp>
#include<multi/adaptors/blas/gemv.hpp>

// namespace multi = boost::multi;

// #include "../../array.hpp"

#include<cmath> // std::isnan
#include<iostream>
#include<algorithm> // std::max

namespace multi = ::boost::multi;
// namespace lapack = ::boost::multi::lapack;

template<class M> decltype(auto) print(M const& C){
	using std::cout;
	using multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout << C[i][j] << ' ';
		cout << std::endl;
	}
	return cout << std::endl;
}

template<class M> decltype(auto) print_1d(M const& C){
	using std::cout;
	using multi::size;
	for(int i = 0; i != size(C); ++i) cout<< C[i] <<' ';
	return cout << std::endl;
}

BOOST_AUTO_TEST_CASE(lapack_geqrf){

	multi::array<double, 2> A = 
		{
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
			{7.0, 8.0, 9.0}
		}
	;
//  multi::lapack::context ctxt;

	multi::array<double, 1> TAU(std::min(size(A), size(~A)));
	multi::array<double, 1> WORK(std::max(1l, 3*size(A)-1));

	multi::lapack::geqrf(ctxt, A, TAU, WORK);

	print(A);
	print(TAU);

}

#if 0
BOOST_AUTO_TEST_CASE(lapack_syev, *boost::unit_test::tolerance(0.00001) ){
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	multi::array<double, 1> WORK(std::max(1l, 3*size(A)-1));
	multi::lapack::syev(multi::blas::filling::upper, A, W, WORK);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::upper, A, W);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::lower, rotated(A), W);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	namespace lapack = multi::lapack;
	multi::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	auto W = lapack::syev(multi::blas::filling::upper, A);
	BOOST_TEST( A[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto A_copy = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( A_copy[2][1] == -0.579092 );
	BOOST_TEST( W[1] == 42.2081 );
}
{
	multi::array<double, 2> A = {
		{167.413, 126.804, 0.},
		{NAN    , 167.381, 0.},
		{NAN    , NAN    , 0.}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( size(A_ref)==3 );
	BOOST_TEST( W[0]==0. );
}
{
	multi::array<double, 2> A = {
		{1. , 1.,  1.},
		{NAN, 2 ,  1.},
		{NAN, NAN, 1.}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
	print(A_ref);
	BOOST_TEST( size(A_ref)==3 );
	BOOST_TEST( W[0]==0. );
}
{
	multi::array<double, 2> A = {{5.}};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[0][0] == 1. );
	BOOST_TEST( W[0]==5. );
}
{
	namespace lapack = multi::lapack;
	multi::array<double, 2> A;
	multi::array<double, 1> W(size(A));
	lapack::syev(lapack::filling::upper, A, W);
}
{
	multi::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto sys = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( sys.eigenvectors[2][1] == -0.579092 );
	BOOST_TEST( sys.eigenvalues[1] == 42.2081 );
}
{
	multi::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto [eigenvecs, eigenvals] = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( eigenvecs[2][1] == -0.579092 );
	BOOST_TEST( eigenvals[1] == 42.2081 );
}
}
#endif

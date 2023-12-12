#ifdef COMPILATION_INSTRUCTIONS
`#nvcc -x cu --expt-relaxed-constexpr`$CXX -D_TEST_MULTI_ADAPTORS_LAPACK_SYEV $0 -o $0x `pkg-config --libs blas lapack` -lboost_unit_test_framework -lcudart -lcusolver&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi lapack adaptor syev"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../lapack/syev.hpp"

#include "../../../array.hpp"
#include "../../../adaptors/cuda.hpp"
#include "../../lapack/cuda.hpp"

BOOST_AUTO_TEST_CASE(multi_lapack_syev, *boost::unit_test::tolerance(0.00001) ){
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
	multi::cuda::managed::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::managed::array<double, 1> W(size(A));
	multi::cuda::managed::array<double, 1> WORK(std::max(1l, 3*size(A)-1));
	multi::lapack::syev(multi::blas::filling::upper, A, W, WORK);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
}
{
	multi::cuda::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::array<double, 1> W(size(A));
	multi::cuda::array<double, 1> WORK(std::max(1l, 3*size(A)-1));
	multi::lapack::syev(multi::blas::filling::upper, A, W, WORK);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
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
	multi::cuda::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::upper, A, W);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
}
{
	multi::cuda::managed::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::managed::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::upper, A, W);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
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
	multi::cuda::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::lower, rotated(A), W);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
}
{
	multi::cuda::managed::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::managed::array<double, 1> W(size(A));
	multi::lapack::syev(multi::blas::filling::lower, rotated(A), W);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
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
	namespace lapack = multi::lapack;
	multi::cuda::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	auto W = lapack::syev(multi::blas::filling::upper, A);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
}
{
	namespace lapack = multi::lapack;
	multi::cuda::managed::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	auto W = lapack::syev(multi::blas::filling::upper, A);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
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
	namespace lapack = multi::lapack;
	multi::cuda::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	auto W = lapack::syev(multi::blas::filling::upper, A);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
}
{
	namespace lapack = multi::lapack;
	multi::cuda::managed::array<double, 2> A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	auto W = lapack::syev(multi::blas::filling::upper, A);
	BOOST_TEST( double(A[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
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
	multi::cuda::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto A_copy = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( double(A_copy[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
}
{
	multi::cuda::managed::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::managed::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto A_copy = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( double(A_copy[2][1]) == -0.579092 );
	BOOST_TEST( double(W[1]) == 42.2081 );
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
	BOOST_TEST( size(A_ref)==3 );
	BOOST_TEST( W[0]==0. );
}
{
	multi::cuda::array<double, 2> A = {
		{1. , 1.,  1.},
		{NAN, 2 ,  1.},
		{NAN, NAN, 1.}
	};
	multi::cuda::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( size(A_ref)==3 );
	BOOST_TEST( double(W[0])==0. );
}
{
	multi::cuda::managed::array<double, 2> A = {
		{1. , 1.,  1.},
		{NAN, 2 ,  1.},
		{NAN, NAN, 1.}
	};
	multi::cuda::managed::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto&& A_ref = lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( size(A_ref)==3 );
	BOOST_TEST( double(W[0])==0. );
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
	multi::cuda::array<double, 2> A = {{5.}};
	multi::cuda::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[0][0] == 1. );
	BOOST_TEST( W[0]==5. );
}
{
	multi::cuda::managed::array<double, 2> A = {{5.}};
	multi::cuda::managed::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	lapack::syev(lapack::filling::upper, A, W);
	BOOST_TEST( A[0][0] == 1. );
	BOOST_TEST( W[0]==5. );
}
{
	multi::array<double, 2> A;
	multi::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	lapack::syev(lapack::filling::upper, A, W);
}
{
	multi::cuda::array<double, 2> A;
	multi::cuda::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	lapack::syev(lapack::filling::upper, A, W);
}
{
	multi::cuda::managed::array<double, 2> A;
	multi::cuda::managed::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
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
	multi::cuda::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto sys = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( double(sys.eigenvectors[2][1]) == -0.579092 );
	BOOST_TEST( double(sys.eigenvalues[1]) == 42.2081 );
}
{
	multi::cuda::managed::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::managed::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto sys = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( double(sys.eigenvectors[2][1]) == -0.579092 );
	BOOST_TEST( double(sys.eigenvalues[1]) == 42.2081 );
}
#if __cpp_structured_bindings
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
{
	multi::cuda::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto [eigenvecs, eigenvals] = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( eigenvecs[2][1] == -0.579092 );
	BOOST_TEST( eigenvals[1] == 42.2081 );
}
{
	multi::cuda::managed::array<double, 2> const A = {
		{167.413, 126.804, 125.114},
		{NAN    , 167.381, 126.746},
		{NAN    , NAN    , 167.231}
	};
	multi::cuda::managed::array<double, 1> W(size(A));
	namespace lapack = multi::lapack;
	auto [eigenvecs, eigenvals] = lapack::syev(lapack::filling::upper, A);
	BOOST_TEST( A[1][2] == 126.746 );
	BOOST_TEST( eigenvecs[2][1] == -0.579092 );
	BOOST_TEST( eigenvals[1] == 42.2081 );
}
#endif
}



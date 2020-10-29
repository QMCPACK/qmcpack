#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&`#nvcc -x cu --expt-relaxed-constexpr`$CXX $0 -o $0x -lcudart -lcublas `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../memory/adaptors/cuda/managed/ptr.hpp"

#include "../../../adaptors/blas.hpp"
#include "../../../adaptors/blas/cuda.hpp"

#include "../../../adaptors/cuda.hpp"
#include "../../../array.hpp"

namespace multi = boost::multi;

template<class M> decltype(auto) print(M const& C){
	using multi::size; using std::cout;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<std::endl;
}

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex, *utf::tolerance(0.00001)){
	using complex = std::complex<double>;
	constexpr complex I{0., 1.};
	multi::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> B = {
		{1. - 9.*I, 3. + 2.*I, 4. + 3.*I},
		{2. - 2.*I, 7. - 2.*I, 1. - 1.*I},
		{3. + 1.*I, 4. + 8.*I, 2. + 7.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_TEST( std::real(B[1][2]) == 2.33846 );
	BOOST_TEST( std::imag(B[1][2]) == -0.0923077 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_rectangular, *utf::tolerance(0.00001)){
	using complex = std::complex<double>;
	constexpr complex I{0., 1.};
	multi::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> B = {
		{1. - 9.*I, 3. + 2.*I},
		{2. - 2.*I, 7. - 2.*I},
		{3. + 1.*I, 4. + 8.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_TEST( real(B[2][0]) == -4.16471 );
	BOOST_TEST( imag(B[2][0]) ==  8.25882 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column, *utf::tolerance(0.00001)){
	using complex = std::complex<double>;
	constexpr complex I{0., 1.};
	multi::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> B = {
		{1. - 9.*I},
		{2. - 2.*I},
		{3. + 1.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;	
	using blas::hermitized;
	trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_TEST( real(B[2][0]) == -4.16471 );
	BOOST_TEST( imag(B[2][0]) ==  8.25882 );
}

using complex = std::complex<double>; complex const I{0, 1};

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column_cpu, *utf::tolerance(0.00001)){
	multi::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::array<complex, 2> const B = {
		{1. - 9.*I},
		{2. - 2.*I},
		{3. + 1.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	auto C = trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_TEST( real(C[2][0]) == -4.16471 );
	BOOST_TEST( imag(C[2][0]) ==  8.25882 );
}

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column_cuda, *utf::tolerance(0.00001)){
	multi::cuda::array<complex, 2> const A = {
		{ 1.,  3.,  4.},
		{NAN,  7.,  1.},
		{NAN, NAN,  8.}
	};
	multi::cuda::array<complex, 2> const B = {
		{1.},
		{2.},
		{3.}
	};
	namespace blas = multi::blas;
	using blas::filling;	
	using blas::hermitized;
	auto Bcpy = trsm(filling::upper, 1., A, B); // B ⬅ α Inv[A].B, B† ⬅ B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	multi::array<complex, 2> Bcpu = Bcpy;
	BOOST_TEST( std::real(Bcpu[2][0]) == 0.375 );
	BOOST_TEST( std::imag(Bcpu[2][0]) == 0.    );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_column_cuda, *utf::tolerance(0.00001)){
	multi::cuda::array<double, 2> const A = {
		{ 1.,  3.,  4.},
		{NAN,  7.,  1.},
		{NAN, NAN,  8.}
	};
	multi::cuda::array<double, 2> B = {
		{1.},
		{2.},
		{3.}
	};
	namespace blas = multi::blas;
	using blas::filling;	
	using blas::hermitized;
	trsm(filling::upper, 1., A, B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	BOOST_TEST( B[2][0] == 0.375 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_column_cuda2, *utf::tolerance(0.00001)){
	multi::cuda::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::cuda::array<complex, 2> B = {
		{1. - 9.*I},
		{2. - 2.*I},
		{3. + 1.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;	
	using blas::hermitized;
	trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	multi::array<complex, 2> Bcpu = B;
	BOOST_TEST( real(Bcpu[2][0]) == -4.16471 );
	BOOST_TEST( imag(Bcpu[2][0]) ==  8.25882 );
}

BOOST_AUTO_TEST_CASE(multi_blas_cuda_trsm_complex, *utf::tolerance(0.00001)){
	multi::cuda::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::cuda::array<complex, 2> const B = {
		{1. - 9.*I, 3. + 2.*I, 4. + 3.*I},
		{2. - 2.*I, 7. - 2.*I, 1. - 1.*I},
		{3. + 1.*I, 4. + 8.*I, 2. + 7.*I}
	};

	namespace blas = multi::blas;
	using blas::filling;	
	using blas::hermitized;
//	auto C = trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
	auto C = trsm(filling::lower, 1., hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit 
}

BOOST_AUTO_TEST_CASE(multi_blas_cuda_managed_trsm_complex, *utf::tolerance(0.00001)){
	multi::cuda::managed::array<complex, 2> const A = {
		{ 1. + 2.*I,  3. - 1.*I,  4. + 9.*I},
		{NAN       ,  7. + 4.*I,  1. + 8.*I},
		{NAN       , NAN       ,  8. + 2.*I}
	};
	multi::cuda::managed::array<complex, 2> const B = {
		{1. - 9.*I, 3. + 2.*I, 4. + 3.*I},
		{2. - 2.*I, 7. - 2.*I, 1. - 1.*I},
		{3. + 1.*I, 4. + 8.*I, 2. + 7.*I}
	};

	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	auto C = trsm(filling::lower, 2.+1.*I, hermitized(A), B); // B=alpha Inv[A†].B, B†=B†.Inv[A], Solve(A†.X=B, X), Solve(X†.A=B†, X), A is upper triangular (with implicit zeros below)
}


#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lcudart -lcublas -lboost_unit_test_framework `pkg-config --libs blas`&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS trsv"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../memory/adaptors/cuda/managed/ptr.hpp"

#include "../../../adaptors/blas/trsv.hpp"
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

using complex = std::complex<double>;
complex const I{0, 1};

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cpu, *utf::tolerance(0.00001)){
	
	multi::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	multi::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(blas::filling::upper, blas::diagonal::general, A, b);
	BOOST_TEST_REQUIRE( real(b[0]) == -1.37259 );
	BOOST_TEST_REQUIRE( real(b[1]) ==  0.2127 );
	BOOST_TEST_REQUIRE( real(b[2]) ==  0.569231 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cuda, *utf::tolerance(0.0001)){
	namespace cuda = multi::cuda;
	cuda::managed::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	cuda::managed::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(blas::filling::upper, blas::diagonal::general, A, b);

	BOOST_TEST_REQUIRE( real(b[0]) == -1.37259  );
	BOOST_TEST_REQUIRE( real(b[1]) ==  0.2127   );
	BOOST_TEST_REQUIRE( real(b[2]) ==  0.569231 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cuda_managed, *utf::tolerance(0.00001)){
	namespace cuda = multi::cuda;
	cuda::managed::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	cuda::managed::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(blas::filling::upper, A, b); // this operation happens in GPU when #include "adaptors/blas/cuda.hpp"

	multi::array<complex, 1> const b_cpu = b;
	BOOST_TEST_REQUIRE( real(b_cpu[0]) == -1.37259  );
	BOOST_TEST_REQUIRE( real(b_cpu[1]) ==  0.2127   );
	BOOST_TEST_REQUIRE( real(b_cpu[2]) ==  0.569231 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_cuda_managed, *utf::tolerance(0.00001)){
	namespace cuda = multi::cuda;
	cuda::managed::array<double, 2> const A = {
		{ 1.,  3.,  4.},
		{NAN,  7.,  1.},
		{NAN, NAN,  8.}
	};
	cuda::managed::array<double, 1> b = {1., 3., 4.};

	blas::trsv(blas::filling::upper, A, b); // this operation happens in GPU when #include "adaptors/blas/cuda.hpp"
	multi::array<double, 1> const b_cpu = b;
	BOOST_TEST_REQUIRE( b_cpu[0] == -2.07143  );
	BOOST_TEST_REQUIRE( b_cpu[1] ==  0.357143  );
	BOOST_TEST_REQUIRE( b_cpu[2] ==  0.5 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cuda2, *utf::tolerance(0.00001)){
	namespace blas = multi::blas;
	multi::cuda::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	multi::cuda::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(blas::filling::upper, blas::diagonal::general, A, b);
	BOOST_TEST_REQUIRE( real(b[0]) == -1.37259 );
	BOOST_TEST_REQUIRE( real(b[1]) ==  0.2127 );
	BOOST_TEST_REQUIRE( real(b[2]) ==  0.569231 );
}


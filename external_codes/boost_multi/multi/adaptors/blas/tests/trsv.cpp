#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&`#nvcc -x cu --expt-relaxed-constexpr`$CXX $0 -o $0x `#-Wno-deprecated-declarations` -lcudart -lcublas -lboost_unit_test_framework `pkg-config --libs blas`&&$0x&&rm $0x $0.cpp; exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS gemm"
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
using blas::filling;
using blas::diagonal;
using blas::transposed;
using blas::hermitized;
using blas::conjugated;
using blas::trsv;

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cpu, *utf::tolerance(0.00001)){
	multi::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	multi::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(filling::upper, diagonal::general, A, b);
	BOOST_TEST( real(b[0]) == -1.37259 );
	BOOST_TEST( real(b[1]) ==  0.2127 );
	BOOST_TEST( real(b[2]) ==  0.569231 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cuda, *utf::tolerance(0.0001)){
	multi::cuda::managed::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	multi::cuda::managed::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(filling::upper, A, b);
	complex b_cpu0 = b[0];
//	std::cout << real(b_cpu0) << std::endl; 
//	BOOST_REQUIRE( std::real(b_cpu0) == -1.37259  );
//	BOOST_REQUIRE( real(b_cpu[1]) ==  0.2127   );
//	BOOST_REQUIRE( real(b_cpu[2]) ==  0.569231 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cuda_managed, *utf::tolerance(0.00001)){
	multi::cuda::managed::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	multi::cuda::managed::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(filling::upper, A, b); // this operation happens in GPU when #include "adaptors/blas/cuda.hpp"
	multi::array<complex, 1> const b_cpu = b;
	std::cout << real(b_cpu[1]) << '\n';
//	BOOST_REQUIRE( real(b_cpu[0]) == -1.37259  );
//	BOOST_REQUIRE( real(b_cpu[1]) ==  0.2127   );
//	BOOST_REQUIRE( real(b_cpu[2]) ==  0.569231 );
}

BOOST_AUTO_TEST_CASE(multi_blas_trsm_double_cuda_managed, *utf::tolerance(0.00001)){
	multi::cuda::managed::array<double, 2> const A = {
		{ 1.,  3.,  4.},
		{NAN,  7.,  1.},
		{NAN, NAN,  8.}
	};
	multi::cuda::managed::array<double, 1> b = {1., 3., 4.};
	blas::trsv(filling::upper, A, b); // this operation happens in GPU when #include "adaptors/blas/cuda.hpp"
	multi::array<double, 1> const b_cpu = b;
	std::cout << b_cpu[1] << '\n';
//	BOOST_REQUIRE( real(b_cpu[0]) == -1.37259  );
//	BOOST_REQUIRE( real(b_cpu[1]) ==  0.2127   );
//	BOOST_REQUIRE( real(b_cpu[2]) ==  0.569231 );
}


#if 0
BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_cuda, *utf::tolerance(0.00001)){
	multi::cuda::array<complex, 2> const A = {
		{ 1. + 1.*I,  3. -  2.*I,  4. + 1.*I},
		{NAN       ,  7. - 10.*I,  1. + 2.*I},
		{NAN       , NAN        ,  8. + 1.*I}
	};
	multi::cuda::array<complex, 1> b = {1. + 2.*I, 3. + 1.*I, 4. + 5.*I};
	blas::trsv(filling::upper, diagonal::general, A, b);
//	BOOST_TEST( real(b[0]) == -1.37259 );
//	BOOST_TEST( real(b[1]) ==  0.2127 );
//	BOOST_TEST( real(b[2]) ==  0.569231 );
}
#endif


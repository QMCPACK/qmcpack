#ifdef COMPILATION_INSTRUCTIONS
nvcc -x cu --expt-relaxed-constexpr`#$CXX` -Wno-deprecated-declarations $0 -o $0x -lcudart -lcublas -lboost_unit_test_framework \
`pkg-config --libs blas` -DBOOST_LOG_DYN_LINK -lboost_system&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS herk"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../adaptors/cuda.hpp" // multi::cuda ns

#include "../../../adaptors/blas/cuda.hpp"
#include "../../../adaptors/blas/herk.hpp"
#include "../../../adaptors/blas/gemm.hpp"

#include "../../../array.hpp"

#include <boost/log/expressions.hpp>

namespace multi = boost::multi;
namespace cuda = multi::cuda;

BOOST_AUTO_TEST_CASE(multi_blas_cuda_herk_complex){
#if 1
	using complex = std::complex<double>;
	complex const I{0, 1};

	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		namespace blas = multi::blas;
		using blas::herk;
		herk(a, c);
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );

		multi::array<complex, 2> const c_copy = herk(1., a);
		BOOST_REQUIRE( c == c_copy );

		using blas::gemm;
		using blas::hermitized;
		BOOST_REQUIRE( gemm(a, hermitized(a)) == herk(a) );
	}
	{
		cuda::array<complex, 2> const acu = a; BOOST_REQUIRE(a == acu);
		cuda::array<complex, 2> ccu({2, 2}, 9999.);
		using multi::blas::herk;
		herk(acu, ccu);
		BOOST_REQUIRE( ccu[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( ccu[0][1] == complex(50., +49.) );

		cuda::array<complex, 2> const ccu_copy = herk(1., acu);
		BOOST_REQUIRE( herk(1., acu) == ccu );
	}
	{
		cuda::managed::array<complex, 2> const amcu = a; BOOST_REQUIRE(a == amcu);
		cuda::managed::array<complex, 2> cmcu({2, 2}, 9999.);
		using multi::blas::herk;
		herk(1., amcu, cmcu);
		BOOST_REQUIRE( cmcu[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( cmcu[0][1] == complex(50., +49.) );

		cuda::managed::array<complex, 2> const cmcu_copy = herk(1., amcu);
		BOOST_REQUIRE( cmcu_copy == cmcu );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		using multi::blas::herk;
		using multi::blas::hermitized;
		herk(1., hermitized(a), c);
		BOOST_REQUIRE( c[2][1] == complex(41, +2) );
		BOOST_REQUIRE( c[1][2] == complex(41, -2) );

		multi::array<complex, 2> const c_copy = herk(1., hermitized(a));
		BOOST_REQUIRE( c_copy == c );
	}
	{
		cuda::array<complex, 2> const acu = a; BOOST_REQUIRE(a == acu);
		cuda::array<complex, 2> ccu({3, 3}, 9999.);
		using multi::blas::herk;
		using multi::blas::hermitized;
		herk(1., hermitized(acu), ccu);
		BOOST_REQUIRE( ccu[2][1] == complex(41, +2) );
		BOOST_REQUIRE( ccu[1][2] == complex(41, -2) );

		cuda::array<complex, 2> const ccu_copy = herk(1., hermitized(acu));
		BOOST_REQUIRE( ccu_copy == ccu );
	}
	{
		cuda::managed::array<complex, 2> const acu = a; BOOST_REQUIRE(a == acu);
		cuda::managed::array<complex, 2> ccu({3, 3}, 9999.);
		using multi::blas::herk;
		using multi::blas::hermitized;
		herk(1., hermitized(acu), ccu);
		BOOST_REQUIRE( ccu[2][1] == complex(41, +2) );
		BOOST_REQUIRE( ccu[1][2] == complex(41, -2) );

		cuda::managed::array<complex, 2> const ccu_copy = herk(1., hermitized(acu));
		BOOST_REQUIRE( ccu_copy == ccu );		
	}
#endif
}

#if 0
BOOST_AUTO_TEST_CASE(multi_blas_cuda_herk_real){
	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	{
		multi::array<double, 2> c({2, 2}, 9999);
		using multi::blas::herk;
		herk(1., a, c);
		BOOST_REQUIRE( c[1][0] == 34 );
		BOOST_REQUIRE( c[0][1] == 34 );

		multi::array<double, 2> const c_copy = herk(1., a);
		BOOST_REQUIRE( c == c_copy );
	}
	{
		cuda::array<double, 2> const acu = a; BOOST_REQUIRE(a == acu);
		cuda::array<double, 2> ccu({2, 2}, 9999.);
		using multi::blas::herk;
		herk(1., acu, ccu);
		BOOST_REQUIRE( ccu[1][0] == 34 );
		BOOST_REQUIRE( ccu[0][1] == 34 );

		cuda::array<double, 2> const ccu_copy = herk(1., acu);
		BOOST_REQUIRE( herk(1., acu) == ccu );
	}
	{
		cuda::array<double, 2> const acu = a; BOOST_REQUIRE(a == acu);
	//	cuda::array<double, 2> ccu({2, 2}, 9999.);
		using multi::blas::herk;
		cuda::array<double, 2> ccu = herk(acu);
		BOOST_REQUIRE( ccu[1][0] == 34 );
		BOOST_REQUIRE( ccu[0][1] == 34 );

		cuda::array<double, 2> const ccu_copy = herk(1., acu);
		BOOST_REQUIRE( herk(1., acu) == ccu );
	}
	{
		cuda::managed::array<double, 2> const amcu = a; BOOST_REQUIRE(a == amcu);
		cuda::managed::array<double, 2> cmcu({2, 2}, 9999.);
		using multi::blas::herk;
		herk(1., amcu, cmcu);
		BOOST_REQUIRE( cmcu[1][0] == 34 );
		BOOST_REQUIRE( cmcu[0][1] == 34 );

		cuda::managed::array<double, 2> const cmcu_copy = herk(1., amcu);
		BOOST_REQUIRE( cmcu_copy == cmcu );
	}
	if(0){
		multi::array<double, 2> c({3, 3}, 9999.);
		using multi::blas::herk;
		using multi::blas::hermitized;
		herk(1., hermitized(a), c);
		BOOST_REQUIRE( c[2][1] == 19 );
		BOOST_REQUIRE( c[1][2] == 19 );

		multi::array<double, 2> const c_copy = herk(1., hermitized(a));
		BOOST_REQUIRE( c_copy == c );
	}
	if(0){
		cuda::array<double, 2> const acu = a; BOOST_REQUIRE(acu == a);
		cuda::array<double, 2> ccu({3, 3}, 9999.);
		using multi::blas::herk;
		using multi::blas::hermitized;
		herk(1., hermitized(acu), ccu);
		BOOST_REQUIRE( ccu[2][1] == 19 );
		BOOST_REQUIRE( ccu[1][2] == 19 );

		cuda::array<double, 2> const c_copy = herk(1., hermitized(a));
		BOOST_REQUIRE( c_copy == ccu );
	}
	if(0){
		cuda::managed::array<double, 2> const amcu = a; BOOST_REQUIRE(amcu == a);
		cuda::managed::array<double, 2> cmcu({3, 3}, 9999.);
		using multi::blas::herk;
		using multi::blas::hermitized;
		herk(1., hermitized(amcu), cmcu);
		BOOST_REQUIRE( cmcu[2][1] == 19 );
		BOOST_REQUIRE( cmcu[1][2] == 19 );

		cuda::managed::array<double, 2> const c_copy = herk(1., hermitized(a));
		BOOST_REQUIRE( c_copy == cmcu );
	}
}
#endif


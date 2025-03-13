#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS herk"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../../adaptors/cuda.hpp" // multi::cuda ns
#include "../../../../adaptors/blas/herk.hpp"


namespace multi = boost::multi;
using complex = std::complex<double>;
complex const I{0, 1};

BOOST_AUTO_TEST_CASE(multi_blas_herk){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	multi::cuda::array<complex, 2> const a_gpu = a;
	namespace blas = multi::blas;
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(1., a, c);
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );

		multi::array<complex, 2> const c_copy = blas::herk(1., a);
		BOOST_REQUIRE( c == c_copy );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		blas::herk(1., blas::H(a), c);
		BOOST_REQUIRE( c[2][1] == complex(41, +2) );
		BOOST_REQUIRE( c[1][2] == complex(41, -2) );

		multi::array<complex, 2> const c_copy = blas::herk(1., blas::H(a));
		BOOST_REQUIRE( c_copy == c );
	}
}


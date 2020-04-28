#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4-*-
clang++ -std=c++14 --cuda-gpu-arch=sm_52 -Ofast -x cuda $0 -o $0x -ltbb -lcudart -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#include "../array.hpp"
#include "../adaptors/cuda.hpp"

#include<complex>
#include<numeric>

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reinterpret array"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

namespace multi = boost::multi;

namespace{

template<class T> struct Complex_{T real; T imag; T const& re() const{return real;}};

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast){
{
	std::complex<double> c{1, 2};
	auto pC = reinterpret_cast<Complex_<double>*>(&c);
	pC->real = 11;
	BOOST_REQUIRE(real(c)==11);
}
{
	multi::array<std::complex<double>, 1> A(10); 
	std::iota( begin(A), end(A), 1.);
	BOOST_REQUIRE( A[8] == 9. );
	auto&& A2 = multi::reinterpret_array_cast<Complex_<double>>(A);
	A2[8].real = 1000.;
	BOOST_REQUIRE( A[8] == 1000. );
}
{
	multi::cuda::managed::array<std::complex<double>, 1> A(10);
	A[8] = std::complex<double>{1000., 2000.};
	auto&& A2 = multi::reinterpret_array_cast<Complex_<double>>(A);
	BOOST_REQUIRE( A[8] == std::complex<double>(1000., 2000.) );

	auto&& AR = multi::member_array_cast<double>(A2, &Complex_<double>::real);
	auto&& AI = multi::member_array_cast<double>(A2, &Complex_<double>::imag);

	BOOST_REQUIRE( AR[8] == 1000. );
	BOOST_REQUIRE( AI[8] == 2000. );

//	auto&& ARP = multi::member_array_cast<double>(A2, &Complex_<double>::re);
}
{
	multi::cuda::array<std::complex<double>, 1> A(10);
	CUDA_SLOW(( A[8] = std::complex<double>{1000., 2000.} ));
	auto&& A2 = multi::reinterpret_array_cast<Complex_<double>>(A);
	BOOST_REQUIRE( CUDA_SLOW( A[8] == std::complex<double>(1000., 2000.) ) );

	auto&& AR = multi::member_array_cast<double>(A2, &Complex_<double>::real);
	auto&& AI = multi::member_array_cast<double>(A2, &Complex_<double>::imag);

	BOOST_REQUIRE( CUDA_SLOW( AR[8] == 1000. ) );
	BOOST_REQUIRE( CUDA_SLOW( AI[8] == 2000. ) );
}
}

}


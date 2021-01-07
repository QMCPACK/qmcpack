#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reinterpret array"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"
//#include "../adaptors/cuda.hpp"

#include<complex>
#include<numeric>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_struct_to_dimension){
	struct vec3{
		double x, y, z;
	};
	multi::array<vec3, 1> A(100);
	A[8] = {1., 2., 3.};
	BOOST_REQUIRE( A[8].y == 2. );

	BOOST_REQUIRE( A.reinterpret_array_cast<double>(3)[8][1] == A[8].y );

	multi::array<double, 2> A2D = A.reinterpret_array_cast<double>(3);
	BOOST_REQUIRE( dimensionality(A2D) == dimensionality(A) + 1 );
	BOOST_REQUIRE( size(A2D) == size(A) );
	BOOST_REQUIRE(  A2D[8][1] ==  A[8].y );
	BOOST_REQUIRE( &A2D[8][1] != &A[8].y );

	A.reinterpret_array_cast<double>(3)[8][1] = 99.;
	BOOST_REQUIRE( A[8].y == 99.);
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_complex_to_real_extra_dimension){
	using complex = std::complex<double>;
	multi::array<complex, 1> A(100, {1, 2});
	BOOST_REQUIRE(  size(A) == 100 );
	BOOST_REQUIRE(( A[0] == complex{1, 2} ));
	
	multi::array<double, 1> B = A.reinterpret_array_cast<double>();
	BOOST_REQUIRE( dimensionality(B) == dimensionality(A) );
	BOOST_REQUIRE( B[0] == 1 and B[1] == 1 );

	multi::array<double, 2> C = A.reinterpret_array_cast<double>(2);	

	BOOST_REQUIRE(( sizes(C)==decltype(sizes(C)){100, 2} ));
	BOOST_REQUIRE( C[5][0] == real(A[5]) );
	BOOST_REQUIRE( C[5][1] == imag(A[5]) );
}

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
	auto&& A2 = A.template reinterpret_array_cast<Complex_<double>>();
	A2[8].real = 1000.;
	BOOST_REQUIRE( A[8] == 1000. );
}

//	BOOST_TEST( size(D) == size(A) );
	
//	auto&& D = A.template reinterpret_array_cast<double[2]>();
//	BOOST_REQUIRE( D[5][0] == real(A[5]) );
//	BOOST_REQUIRE( D[5][1] == imag(A[5]) );

//	multi::array<double, 2> C = A.template reinterpret_array_cast<double>(2);
//	multi::array<double, 2> C = reinterpret_array_cast<double>(A, 2);
	
//	multi::array<double, 2> C = A.template reinterpret_array_cast<double, 2>();
//	multi::array<double, 2> C = reinterpret_array_cast<double>(A, 2);

}

#if 0
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
#endif
#if 0
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
#endif


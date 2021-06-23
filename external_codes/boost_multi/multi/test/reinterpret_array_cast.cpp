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

	BOOST_REQUIRE( & A[8].x == & A.reinterpret_array_cast<double>(3)[8][0] );
	BOOST_REQUIRE( & A[8].y == & A.reinterpret_array_cast<double>(3)[8][1] );
	BOOST_REQUIRE( & A[8].z == & A.reinterpret_array_cast<double>(3)[8][2] );
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_complex_to_real_extra_dimension){
	using complex = std::complex<double>;
	multi::array<complex, 1> A(100, complex{1, 2});
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

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_tuple_as_extra_dimension){
	using vector3 = std::array<double, 3>;
//	using vector3 = std::tuple<double, double, double>; // for tuples reinterpret_array_cast is implementation dependent!!

	vector3 v;
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast, cppcoreguidelines-avoid-c-arrays, hicpp-avoid-c-arrays, modernize-avoid-c-arrays): test
	BOOST_REQUIRE( &reinterpret_cast<double(&)[3]>(v)[1] == &std::get<1>(v) ); 
	{
		multi::array<vector3, 1> A(10);
		BOOST_REQUIRE( &A.reinterpret_array_cast<double>(3)[2][1] == &std::get<1>(A[2]) );
	}
	{
		multi::array<vector3, 2> A({10, 20});
		BOOST_REQUIRE( &A.reinterpret_array_cast<double>(3)[5][7][2] == &std::get<2>(A[5][7]) );
	}

	{
		multi::array<vector3, 2> const A({4, 5}, vector3{1., 2., 3.});

		BOOST_REQUIRE( dimensionality(A.reinterpret_array_cast<double>(3)) == 3 );
		BOOST_REQUIRE( A.reinterpret_array_cast<double>(3).num_elements() == A.num_elements()*3 );
		BOOST_REQUIRE( A.reinterpret_array_cast<double>(3).size() == 4 );
		BOOST_REQUIRE( A.reinterpret_array_cast<double>(3)[0].size() == 5 );
		BOOST_REQUIRE( A.reinterpret_array_cast<double>(3)[0][0].size() == 3 );
		BOOST_REQUIRE( &A.reinterpret_array_cast<double>(3)[2][3][0] == &std::get<0>(A[2][3]) );
		BOOST_REQUIRE( &A.reinterpret_array_cast<double>(3)[2][3][1] == &std::get<1>(A[2][3]) );
		BOOST_REQUIRE( &A.reinterpret_array_cast<double>(3)[2][3][2] == &std::get<2>(A[2][3]) );

		multi::array<double, 3> const B = A.reinterpret_array_cast<double>(3);
		BOOST_REQUIRE( B[2][3][0] == std::get<0>(A[2][3]) );
		BOOST_REQUIRE( B[2][3][1] == std::get<1>(A[2][3]) );
		BOOST_REQUIRE( B[2][3][2] == std::get<2>(A[2][3]) );

		auto C = +A.reinterpret_array_cast<double>(3);
		BOOST_REQUIRE( C == B );
	}
}

template<class T> struct Complex_{T real; T imag;};

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast){
{
	std::complex<double> c{1, 2};
	auto *pC = reinterpret_cast<Complex_<double>*>(&c); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	pC->real = 11;
	BOOST_REQUIRE(real(c)==11);
}
{
	multi::array<std::complex<double>, 1> A(10); 
	std::iota( begin(A), end(A), 1.);
	BOOST_REQUIRE( A[8] == 9. );
	auto&& A2 = A.reinterpret_array_cast<Complex_<double>>();
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

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_realcomplex){
	using complex = std::complex<double>;
{
	complex c{1, 2};
	auto *pC = reinterpret_cast<std::array<double, 2>*>(&c); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	(*pC)[0] = 11;
	BOOST_REQUIRE( pC );
	BOOST_REQUIRE(real(c)==11);
}
{
	complex c{1, 2};
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast, cppcoreguidelines-avoid-c-arrays, hicpp-avoid-c-arrays, modernize-avoid-c-arrays): test purposes
	auto *pC = reinterpret_cast<double(*)[2]>(&c); 
	(*pC)[0] = 11;
	BOOST_REQUIRE( pC );
	BOOST_REQUIRE(real(c)==11);
}
{
	multi::array<complex, 1> A(10); 
	auto&& A2 = A.reinterpret_array_cast<double>(2);
	A2[8][0] = 1000.;
	A2[8][1] = 2000.;
	BOOST_REQUIRE( A[8] == std::complex<double>(1000., 2000.) );
}
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_pair_to_complex){
	using complex = std::complex<double>;
	using pair    = std::pair<double, double>;
	multi::array<complex, 2> A({10, 10}, complex{3., 4.});

	multi::array<complex, 2> const& Aconst = A;
	auto&& A_block = Aconst({0, 5}, {0, 5});

	auto&& Apair_block = A_block.template reinterpret_array_cast<pair const>(); // const is important // cppcheck needs `template`
	
	BOOST_REQUIRE( (void const*)&Apair_block[1][2] == (void*)&A[1][2] );
	
	auto&& Adoubles_block = A_block.reinterpret_array_cast<double const>(2);
	BOOST_REQUIRE( (void const*)&Adoubles_block[1][2][0] == (void*)&A[1][2] );
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


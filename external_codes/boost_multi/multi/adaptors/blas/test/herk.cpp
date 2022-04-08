#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lcudart -lcublas -lboost_unit_test_framework `pkg-config --libs blas`&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS herk"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

//#include "../../../adaptors/cuda.hpp" // multi::cuda ns

//#include "../../../adaptors/blas/cuda.hpp"
#include "../../../adaptors/blas/gemm.hpp"
#include "../../../adaptors/blas/herk.hpp"
#include "../../../adaptors/blas/nrm2.hpp"

#include "../../../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_herk) {
	namespace blas = multi::blas;
	using complex = std::complex<double>; constexpr complex I{0, 1};

	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};

	 {
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(a, c);
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );

		multi::array<complex, 2> const c_copy = blas::herk(1., a);
		BOOST_REQUIRE( c == c_copy );

		BOOST_REQUIRE( +blas::gemm(1., a, blas::H(a)) == blas::herk(a) );
	}
}

BOOST_AUTO_TEST_CASE(inq_case){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{0.,  1.,  2.},
		{3.,  4.,  5.},
		{6.,  7.,  8.},
		{9., 10., 11.}
	};
	BOOST_REQUIRE( (+blas::gemm(1., a, blas::T(a)))[1][2] == 86. );
	{
		multi::array<double, 2> c({4, 4});
		blas::herk(1.0, a, c);
		BOOST_REQUIRE( c[1][2] == (+blas::gemm(1., a, blas::T(a)))[1][2] );
	//  BOOST_REQUIRE( c[2][1] == (+blas::gemm(1., a, blas::T(a)))[2][1] );
	}
	{
		multi::array<double, 2> c = blas::herk(1.0, a);
		BOOST_REQUIRE( c == +blas::gemm(1., a, blas::T(a)) );
	}
	{
		BOOST_REQUIRE( blas::herk(a) == +blas::gemm(1., a, blas::T(a)) );
	}
	{
		BOOST_REQUIRE( blas::herk(2.0, a) == +blas::gemm(2.0, a, blas::T(a)) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_real){
	namespace blas = multi::blas;
	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	{
		multi::array<double, 2> c({2, 2}, 9999);
		blas::herk(1., a, c);
//		BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case){
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{1., 2., 3.}};
	multi::array<double, 2> B = blas::herk(A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_case_scale){
	namespace blas = multi::blas;
	multi::array<double, 2> const A = {{1., 2., 3.}};
	multi::array<double, 2> B = blas::herk(0.1, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_TEST( B[0][0] == (1.*1. + 2.*2. + 3.*3.)*0.1 );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case){
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const A = { {1., 2., 3.} };
	multi::array<complex, 2> B = blas::herk(1.0, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_real_case_scale, *boost::unit_test::tolerance(0.00001)){
	namespace blas = multi::blas;
	using complex = std::complex<double>;
	multi::array<complex, 2> const A = {{1., 2., 3.}};
	multi::array<complex, 2> B = blas::herk(0.1, A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_TEST( real( B[0][0]/0.1 ) == 1.*1. + 2.*2. + 3.*3. );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case){
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> const A = {{1. + 2.*I, 2.+3.*I, 3. + 4.*I}};
	multi::array<complex, 2> B = blas::herk(A);
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(A)[0][0])) == blas::nrm2(A[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_out_param){
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> const A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	multi::array<complex, 2> B({1, 1});
	BOOST_REQUIRE( size(B) == 1 );

	blas::herk(blas::filling::upper, 1.0, blas::H(A), 0.0, B);

	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(B[0][0])) == blas::nrm2(blas::T(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized){
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};

	namespace blas = multi::blas;
	multi::array<complex, 2> B = blas::herk(blas::H(A));

	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(A))[0][0])) == blas::nrm2(rotated(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk1x1_complex_case_hermitized_auto){
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> A = {{1. + 2.*I}, {2.+3.*I}, {3. + 4.*I}};
	auto B = blas::herk(1., blas::hermitized(A));
	static_assert( std::is_same<decltype(B), multi::array<complex, 2>>{}, "!" );
	BOOST_REQUIRE( size(B) == 1 );
	BOOST_REQUIRE( B[0][0] == std::norm(1. + 2.*I) + std::norm(2.+3.*I) + std::norm(3. + 4.*I) );

	BOOST_TEST( std::sqrt(real(blas::herk(blas::H(A))[0][0])) == blas::nrm2(rotated(A)[0])() );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_identity){
	namespace blas = multi::blas;
	using complex = std::complex<double>; auto const I = complex{0., 1.};

	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};

	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(blas::filling::lower, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		static_assert(blas::is_conjugated<decltype(blas::H(c))>{}, "!" );

		blas::herk(blas::filling::lower, 1., a, 0., blas::H(c)); // c†=c=aa†=(aa†)†, `c` in upper triangular

		BOOST_REQUIRE( blas::H(c)[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( blas::H(c)[0][1]==9999. );
	}
	{
	//	multi::array<complex, 2> c({2, 2}, 9999.);
	//	blas::herk(blas::filling::lower, 1., a, 0., blas::T(c)); // c†=c=aa†=(aa†)†, `c` in lower triangular
	//	BOOST_REQUIRE( transposed(c)[1][0]==complex(50., -49.) );
	//	BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., transposed(a), 0., c); // c†=c=aT(aT)† not supported
	//	print(c);
	//	BOOST_REQUIRE( c[1][0]==complex(52., -90.) );
	//	BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., transposed(a), 0., hermitized(c)); // c†=c=aT(aT)† not supported
	//	BOOST_REQUIRE( hermitized(c)[1][0]==complex(52., -90.) );
	//	BOOST_REQUIRE( hermitized(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(blas::filling::lower, 1., blas::T(a), 0., blas::T(c)); // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( transposed(c)[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		blas::herk(blas::filling::lower, 1., blas::T(a), 0., blas::H(blas::T(c))); // c†=c=aT(aT)† not supported
		BOOST_REQUIRE( blas::H(blas::T(c))[1][0]==complex(52., -90.) );
		BOOST_REQUIRE( blas::H(blas::T(c))[0][1]==9999. );
	}
	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		using namespace multi::blas;
//		blas::herk(blas::filling::lower, 1., blas::T(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
//		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(blas::U, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		blas::herk(1., a, c); // c†=c=aa†=(aa†)†
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		blas::herk(blas::L, 1., blas::H(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(52., 90.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
	//	multi::array<complex, 2> c({3, 3}, 9999.);
	//	using namespace multi::blas;
	//	herk(filling::lower, 1., transposed(a), 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
	//	BOOST_REQUIRE( c[0][1]==9999. );
	//	BOOST_REQUIRE( c[1][0]==complex(52., 90.) );
	}
}

#if 0

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_real_case){
	multi::array<complex, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::transposed;
	using blas::hermitized;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);

		herk(filling::lower, 1., hermitized(a), 0., c);//c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(19.,0.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), 0., c);//c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][2]==complex(19.,0.) );
		BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::upper, 1., hermitized(a), 0., transposed(c));//c†=c=a†a=(a†a)†, `c` in lower triangular
	//	print(transposed(c));
	//	BOOST_REQUIRE( c[1][2]==complex(19.,0.) );
	//	BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		using blas::transposed;
	//	herk(filling::upper, 1., transposed(a), 0., c);//c_†=c_=a_†a_=(a_†a_)†, `c_` in lower triangular
	//	BOOST_REQUIRE( c[2][1] == 9999. );
	//	BOOST_REQUIRE( c[1][2] == 19. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_basic_transparent_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, information in `c` lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		using multi::blas::herk;
		herk(filling::upper, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, `c` in upper triangular
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
		BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		herk(filling::lower, 1., a, 0., c); // c†=c=aa†, `a` and `c` are c-ordering, information in c lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		herk(filling::upper, 1., a, 0., c); //c†=c=aa†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_basic_enum_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::hermitized;
	using blas::transposed;
	{
	//	multi::array<complex, 2> c({2, 2}, 8888.);
	//	std::cerr << "here" << std::endl;
	//	herk(filling::lower, 1., hermitized(transposed(a)), 0., c); //c†=c=a†a=(a†a)†, `c` in lower triangular
	//	print(c) << std::endl;
	//	std::cerr << "there" << std::endl;
	//	BOOST_REQUIRE( c[0][1]==complex(41.,2.) );
	//	BOOST_REQUIRE( c[1][0]==8888. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c); //c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		using namespace multi::blas;
		herk(filling::upper, 1., hermitized(a), 0., c); //c†=c=a†a=(a†a)†, `c` in upper triangular
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
		BOOST_REQUIRE( c[2][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using namespace multi::blas;
		herk(filling::lower, 1., a, 0., c); // c†=c=aa†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using namespace multi::blas;
		herk(filling::upper, 1., a, 0., c); // c†=c=aa†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_basic_explicit_enum_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	BOOST_REQUIRE( herk(hermitized(a)) == gemm(hermitized(a), a) );
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
	//	herk(filling::lower, 1., hermitized(a), 0., transposed(c)); // c†=c=a†a=(a†a)†, `c` in lower triangular
	//	print(transposed(c));
	//	BOOST_REQUIRE( c[2][1]==complex(41.,2.) );
	//	BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::lower, 1., hermitized(transposed(a)), 0., transposed(c)); // c†=c=a†a=(a†a)†, `c` in lower triangular
		BOOST_REQUIRE( transposed(c)[1][0]==complex(50.,+49.) );
		BOOST_REQUIRE( transposed(c)[0][1]==9999. );
	}
//	BOOST_REQUIRE( herk(hermitized(transposed(a))) == gemm(hermitized(transposed(a)), transposed(a)) );
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), 0., c); // c†=c=a†a=(a†a)†, `c` in upper triangular
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
		BOOST_REQUIRE( c[2][1]==9999. );
		BOOST_REQUIRE( herk(hermitized(a)) == gemm(hermitized(a), a) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::lower, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
		BOOST_REQUIRE( herk(a) == gemm(a, hermitized(a)) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
		BOOST_REQUIRE( herk(a) == gemm(a, hermitized(a)) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 2., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(100., 98.) );
		BOOST_REQUIRE( c[1][0]==9999. );
		BOOST_REQUIRE( herk(2., a) == gemm(2., a, hermitized(a)) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, 0., c); // c†=c=aa†=(aa†)†, `c` in upper triangular
		BOOST_REQUIRE( c[0][1]==complex(50., 49.) );
		BOOST_REQUIRE( c[1][0]==9999. );
		BOOST_REQUIRE( herk(1., a) == gemm(1., a, hermitized(a)) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_operator_interface){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		namespace blas = multi::blas;
		using blas::filling;
		using blas::hermitized;
		herk(filling::lower, 1., hermitized(a), 0., c); // c=c†=a†a, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41., 2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi:: blas::filling;
		herk(filling::lower, 1., a, 0., c); // c=c†=aa†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		herk(1., a, c); // c=c†=aa†
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==complex(50., +49.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		namespace blas = multi::blas;
		using blas::filling;
		using blas::hermitized;
		herk(filling::lower, 1., hermitized(a), 0., c); // c=c†=a†a, `c` in lower triangular
		herk(filling::upper, 1., hermitized(a), 0., c);
		BOOST_REQUIRE( c[2][1]==complex(41., 2.) );
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_operator_interface_implicit_no_sum){
	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		namespace blas = multi::blas;
		using blas::filling;
		using blas::hermitized;
		herk(filling::lower, 1., hermitized(a), c); // c=c†=a†a, `c` in lower triangular
		BOOST_REQUIRE( c[2][1]==complex(41., 2.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		using multi::blas::filling;
		herk(filling::lower, 1., a, c); // c=c†=aa†, `c` in lower triangular
		BOOST_REQUIRE( c[1][0]==complex(50., -49.) );
		BOOST_REQUIRE( c[0][1]==9999. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_ordering_and_symmetrization){

	multi::array<complex, 2> const a = {
		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
	};
	namespace blas = multi::blas;
	using blas::herk;
	using blas::hermitized;
	using blas::filling;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
		BOOST_REQUIRE( c[2][1]==9999. );
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(1., hermitized(a), c); // c†=c=a†a
		BOOST_REQUIRE( c[2][1]==complex(41., +2.) );
		BOOST_REQUIRE( c[1][2]==complex(41., -2.) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa† // c implicit hermitic in upper
		BOOST_REQUIRE( c[1][0] == 9999. );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(1., a, c); // c†=c=aa†
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, 1., a); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(1., a); // c†=c=aa†
		BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::hermitized;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(filling::upper, 1., hermitized(a)); // c†=c=a†a

		BOOST_REQUIRE( size(hermitized(a))==3 );
	//	BOOST_REQUIRE( c[2][1] == complex(41., +2.) );
		BOOST_REQUIRE( c[1][2] == complex(41., -2.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(filling::upper, a); // c†=c=a†a
//		what(multi::pointer_traits<decltype(base(a))>::default_allocator_of(base(a)));
	//	BOOST_REQUIRE( c[1][0] == complex(50., -49.) );
		BOOST_REQUIRE( c[0][1] == complex(50., +49.) );
	}
	{
		using multi::blas::herk;
		using multi::blas::hermitized;
		using multi::blas::filling;
		multi::array<complex, 2> c = herk(filling::upper, hermitized(a)); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1] == complex(41., +2.) );
		BOOST_REQUIRE( c[1][2] == complex(41., -2.) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_size1_real_case){
	multi::array<complex, 2> const a = {
		{1., 3., 4.}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({1, 1}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa†
		BOOST_TEST( c[0][0] == 26. );
	}
	BOOST_TEST( herk(a) == gemm(a, hermitized(a)) );
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_size1){
	multi::array<complex, 2> const a = {
		{1. + 4.*I, 3. + 2.*I, 4. - 1.*I}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({1, 1}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa†
		BOOST_TEST( c[0][0] == 47. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_size0){
	multi::array<complex, 2> const a;
	using namespace multi::blas;
	{
		multi::array<complex, 2> c;
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_TEST( c[0][0] == 47. );
	}
}


BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_automatic_ordering_and_symmetrization_real_case){

	multi::array<complex, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	using namespace multi::blas;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
	{
		multi::array<complex, 2> c({2, 2}, 9999.);
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, 1., a); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, 1., hermitized(a)); // c†=c=a†a
		BOOST_REQUIRE( size(hermitized(a))==3 );
	//	BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, a); // c†=c=a†a
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<complex, 2> c = herk(filling::upper, hermitized(a)); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
}


BOOST_AUTO_TEST_CASE(multi_blas_herk_real_automatic_ordering_and_symmetrization_real_case){

	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	{
		multi::array<double, 2> c({3, 3}, 9999.);
		using multi::blas::hermitized;
		using multi::blas::herk;
		using multi::blas::filling;
	//	herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
	//	BOOST_REQUIRE( c[2][1]==19. );
	//	BOOST_REQUIRE( c[1][2]==19. );
	}
	{
		multi::array<double, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		using multi::blas::filling;
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		multi::array<double, 2> c({2, 2}, 9999.);
		using multi::blas::herk;
		using multi::blas::filling;
		herk(filling::upper, 1., a, c); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		using multi::blas::herk;
		using multi::blas::filling;
		multi::array<double, 2> c = herk(filling::upper, 1., a); // c†=c=aa†
	//	BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		using multi::blas::herk;
		multi::array<complex, 2> c = herk(a); // c†=c=a†a
		BOOST_REQUIRE( c[1][0] == 34. );
		BOOST_REQUIRE( c[0][1] == 34. );
	}
	{
		using multi::blas::herk;
		using multi::blas::hermitized;
		multi::array<complex, 2> c = herk(hermitized(a)); // c†=c=a†a
		BOOST_REQUIRE( c[2][1]==19. );
		BOOST_REQUIRE( c[1][2]==19. );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_real_case){
	multi::array<double, 2> const a = {
		{ 1., 3., 4.},
		{ 9., 7., 1.}
	};
	using multi::blas::filling;
	{
		static_assert( not boost::multi::blas::is_complex_array<multi::array<double, 2>>{} , "!");
		multi::array<double, 2> c({2, 2}, 9999.);
		syrk(filling::lower, 1., a, 0., c);//c†=c=aa†=(aa†)†, `c` in lower triangular
	}
	{
		multi::array<double, 2> c({2, 2}, 9999.);
		herk(filling::lower, 1., a, 0., c);//c†=c=aa†=(aa†)†, `c` in lower triangular
	}
	{
		static_assert( not boost::multi::blas::is_complex_array<multi::array<double, 2>>{} , "!");
		multi::array<double, 2> c = herk(filling::upper, a);//c†=c=aa†=(aa†)†, `c` in lower triangular
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_real_case_1d){
	multi::array<complex, 2> const a = {
		{ 1., 3., 4.},
	};
	namespace blas = multi::blas;
	using blas::filling;
	using blas::transposed;
	using blas::hermitized;
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(filling::lower, 1., hermitized(a), 0., c);//c†=c=a†a=(a†a)†, `c` in lower triangular
		print(c);
		BOOST_REQUIRE( c[2][1]==complex(12.,0.) );
		BOOST_REQUIRE( c[1][2]==9999. );
	}
	{
		multi::array<complex, 2> c({3, 3}, 9999.);
		herk(2., hermitized(a), c);//c†=c=a†a=(a†a)†, `c` in lower triangular

		BOOST_REQUIRE( c[2][1]==complex(24.,0.) );
		BOOST_REQUIRE( c[1][2]==complex(24.,0.) );
		multi::array<complex, 2> c_gemm({3, 3});
	//	gemm(2., hermitized(a), a, c_gemm);
	}
}
#endif

#if 0
BOOST_AUTO_TEST_CASE(multi_blas_herk_complex_timing){
	multi::array<complex, 2> const a({4000, 4000}); std::iota(data_elements(a), data_elements(a) + num_elements(a), 0.2);
	multi::array<complex, 2> c({4000, 4000}, 9999.);
	boost::timer::auto_cpu_timer t;
	using multi::blas::herk;
	using multi::blas::hermitized;
	using multi::blas::filling;
	herk(filling::upper, 1., hermitized(a), c); // c†=c=a†a
}
#endif


//BOOST_AUTO_TEST_CASE(multi_blas_cuda_herk_complex){
//	namespace blas = multi::blas;
//	multi::array<complex, 2> const a = {
//		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
//		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
//	};
//	{
//		cuda::array<complex, 2> const acu = a;
//		BOOST_REQUIRE(a == acu);

//		cuda::array<complex, 2> ccu({2, 2}, 9999.);
//		blas::herk(acu, ccu);
//		BOOST_REQUIRE( ccu[1][0] == complex(50., -49.) );
//		BOOST_REQUIRE( ccu[0][1] == complex(50., +49.) );

//		cuda::array<complex, 2> const ccu_copy = blas::herk(1., acu);
//		BOOST_REQUIRE( blas::herk(1., acu) == ccu );
//	}
//	{
//		cuda::managed::array<complex, 2> const amcu = a; BOOST_REQUIRE(a == amcu);
//		cuda::managed::array<complex, 2> cmcu({2, 2}, 9999.);

//		blas::herk(1., amcu, cmcu);
//		BOOST_REQUIRE( cmcu[1][0] == complex(50., -49.) );
//		BOOST_REQUIRE( cmcu[0][1] == complex(50., +49.) );

//		cuda::managed::array<complex, 2> const cmcu_copy = blas::herk(1., amcu);
//		BOOST_REQUIRE( cmcu_copy == cmcu );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		blas::herk(1., blas::H(a), c);
//		BOOST_REQUIRE( c[2][1] == complex(41, +2) );
//		BOOST_REQUIRE( c[1][2] == complex(41, -2) );

//		multi::array<complex, 2> const c_copy = blas::herk(1., blas::H(a));
//		BOOST_REQUIRE( c_copy == c );
//	}
//	{
//		cuda::array<complex, 2> const acu = a;
//		BOOST_REQUIRE(a == acu);

//		cuda::array<complex, 2> ccu({3, 3}, 9999.);

//		blas::herk(1., blas::H(acu), ccu);
//		BOOST_REQUIRE( ccu[2][1] == complex(41, +2) );
//		BOOST_REQUIRE( ccu[1][2] == complex(41, -2) );

//		cuda::array<complex, 2> const ccu_copy = blas::herk(1., blas::H(acu));
//		BOOST_REQUIRE( ccu_copy == ccu );
//	}
//	{
//		cuda::managed::array<complex, 2> const acu = a; BOOST_REQUIRE(a == acu);
//		cuda::managed::array<complex, 2> ccu({3, 3}, 9999.);

//		blas::herk(1., blas::H(acu), ccu);
//		BOOST_REQUIRE( ccu[2][1] == complex(41, +2) );
//		BOOST_REQUIRE( ccu[1][2] == complex(41, -2) );

//		cuda::managed::array<complex, 2> const ccu_copy = blas::herk(1., blas::H(acu));
//		BOOST_REQUIRE( ccu_copy == ccu );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_cuda_herk_n_complex){
//	namespace blas = multi::blas;
//	multi::array<complex, 2> const a = {
//		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
//		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
//	};
//	blas::context ctxt;
//	{
//		multi::array<complex, 2> c({2, 2}, 9999.);
//		blas::herk_n(ctxt, blas::filling::upper, 1., a.begin(), a.size(), 0., c.begin());
//		BOOST_TEST_REQUIRE( c[0][1] == complex(50., +49.) );
//		BOOST_TEST_REQUIRE( c[1][0] == 9999. );
//	}
//	{
//		multi::array<complex, 2> c({2, 2}, 9999.);
//		blas::herk_n(ctxt, blas::filling::lower, 1., a.begin(), a.size(), 0., c.begin());
//		BOOST_TEST_REQUIRE( c[0][1] == 9999. );
//		BOOST_TEST_REQUIRE( c[1][0] == complex(50., -49.) );
//	}
//	{
//		multi::array<complex, 2> c({2, 2}, 9999.);
//		blas::herk_n(ctxt, blas::filling::lower, 1., a.begin(), a.size(), 0., c.begin());
//		blas::herk_n(ctxt, blas::filling::upper, 1., a.begin(), a.size(), 0., c.begin());
//		BOOST_TEST_REQUIRE( c[0][1] == complex(50., +49.) );
//		BOOST_TEST_REQUIRE( c[1][0] == complex(50., -49.) );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		blas::herk_n(ctxt, blas::filling::lower, 1., blas::H(a).begin(), blas::H(a).size(), 0., c.begin());
//		BOOST_TEST_REQUIRE( c[1][2] == 9999. );
//		BOOST_TEST_REQUIRE( c[2][1] == complex(41., +2.) );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		blas::herk_n(ctxt, blas::filling::upper, 1., blas::H(a).begin(), blas::H(a).size(), 0., c.begin());
//		BOOST_TEST_REQUIRE( c[1][2] == complex(41., -2.) );
//		BOOST_TEST_REQUIRE( c[2][1] == 9999. );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		blas::herk_n(ctxt, blas::filling::lower, 1., blas::H(a).begin(), blas::H(a).size(), 0., c.begin());
//		blas::herk_n(ctxt, blas::filling::upper, 1., blas::H(a).begin(), blas::H(a).size(), 0., c.begin());
//		BOOST_TEST_REQUIRE( c[1][2] == complex(41., -2.) );
//		BOOST_TEST_REQUIRE( c[2][1] == complex(41., +2.) );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		blas::herk_n(ctxt, 1., blas::H(a).begin(), blas::H(a).size(), c.begin());
//		BOOST_TEST_REQUIRE( c[1][2] == complex(41., -2.) );
//		BOOST_TEST_REQUIRE( c[2][1] == complex(41., +2.) );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_cuda_herk_row){
//	namespace blas = multi::blas;
//	auto const a = []{
//		multi::array<complex, 2> ret({1, 100});
//		std::generate(begin(ret[0]), end(ret[0]), [c=complex{1, 2}]()mutable{return c+=2.;});
//		return ret;
//	}();
//	BOOST_REQUIRE( size(a) == 1 );
//	{
//		BOOST_REQUIRE( +blas::gemm(1., a, blas::H(a)) == blas::herk(a) );

//		cuda::array<complex, 2> const agpu = a;
//		BOOST_REQUIRE( blas::gemm(agpu, blas::H(agpu)) == blas::herk(agpu) );

//		cuda::managed::array<complex, 2> const amng = a;
//		BOOST_REQUIRE( blas::gemm(amng, blas::H(amng)) == blas::herk(amng) );
//	}
//}

//#if 1
//BOOST_AUTO_TEST_CASE(multi_blas_cuda_herk_real){
//	namespace blas = multi::blas;
//	multi::array<double, 2> const a = {
//		{ 1., 3., 4.},
//		{ 9., 7., 1.}
//	};
//	{
//		multi::array<double, 2> c({2, 2}, 9999);
//		blas::herk(1., a, c);
//		BOOST_REQUIRE( c[1][0] == 34 );
//		BOOST_REQUIRE( c[0][1] == 34 );

//	//	multi::array<double, 2> const c_copy = blas::herk(1., a);
//	//	BOOST_REQUIRE( c == c_copy );
//	}
//	{
//		cuda::array<double, 2> acu = a;
//		BOOST_REQUIRE(a == acu);

//		cuda::array<double, 2> ccu({2, 2}, 9999.);

//	//	blas::herk(acu, ccu);
//	//	BOOST_REQUIRE( ccu[1][0] == 34 );
//	//	BOOST_REQUIRE( ccu[0][1] == 34 );

//	//	cuda::array<double, 2> const ccu_copy = blas::herk(1., acu);
//	//	BOOST_REQUIRE( herk(1., acu) == ccu );
//	}

//}
//#endif

#if 0
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
	if(0) {
		multi::array<double, 2> c({3, 3}, 9999.);
		using multi::blas::herk;
		using multi::blas::hermitized;
		herk(1., hermitized(a), c);
		BOOST_REQUIRE( c[2][1] == 19 );
		BOOST_REQUIRE( c[1][2] == 19 );

		multi::array<double, 2> const c_copy = herk(1., hermitized(a));
		BOOST_REQUIRE( c_copy == c );
	}
	if(0) {
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
	if(0) {
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


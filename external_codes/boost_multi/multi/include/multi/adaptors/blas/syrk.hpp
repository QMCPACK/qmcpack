#ifndef MULTI_ADAPTORS_BLAS_SYRK_HPP  // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
#define MULTI_ADAPTORS_BLAS_SYRK_HPP
// Copyright 2019-2021 Alfredo A. Correa

#include "../blas/core.hpp"
#include "../blas/filling.hpp"
#include "../blas/numeric.hpp"

namespace boost::multi::blas {

using core::syrk;

template<typename AA, typename BB, class A2D, class C2D>
auto syrk(filling c_side, AA alpha, A2D const& a, BB beta, C2D&& c) {  // NOLINT(readability-identifier-length) BLAS naming
//->decltype(syrk('\0', '\0', size(c), size(a), alpha, base(a), stride(rotated(a)), beta, base(c), stride(c)), std::forward<C2D>(c)){
	assert( size(c) == size(rotated(c)) );
	if(stride(a)==1) {
		if(stride(c)==1) {syrk(flip(c_side)==filling::upper?'L':'U', 'N', size(c), size(a         ), &alpha, base(a), stride(rotated(a)), &beta, base(c), stride(rotated(c)));}
		else             {syrk(c_side      ==filling::upper?'L':'U', 'N', size(c), size(rotated(a)), &alpha, base(a), stride(rotated(a)), &beta, base(c), stride(        c ));}
	} else {
		if(stride(c)==1) {syrk(flip(c_side)==filling::upper?'L':'U', 'T', size(c), size(rotated(a)), &alpha, base(a), stride(a), &beta, base(c), stride(rotated(c)));}
		else             {syrk(c_side      ==filling::upper?'L':'U', 'T', size(c), size(rotated(a)), &alpha, base(a), stride(a), &beta, base(c), stride(        c ));}
	}
	return std::forward<C2D>(c);
}

template<typename AA, class A2D, class C2D>
auto syrk(filling c_side, AA alpha, A2D const& a, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(syrk(c_side, alpha, a, 0., std::forward<C2D>(c))) {
	return syrk(c_side, alpha, a, 0., std::forward<C2D>(c)); }

template<typename AA, class A2D, class C2D>
auto syrk(AA alpha, A2D const& a, C2D&& c)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(syrk(filling::upper, alpha, a, syrk(filling::lower, alpha, a, std::forward<C2D>(c)))) {
	return syrk(filling::upper, alpha, a, syrk(filling::lower, alpha, a, std::forward<C2D>(c))); }

template<typename AA, class A2D, class Ret = typename A2D::decay_type>
[[nodiscard]]  // ("because input argument is const")
// this decay in the return type is important
// NOLINTNEXTLINE(readability-identifier-length) BLAS naming
auto syrk(AA alpha, A2D const& a) -> std::decay_\
t<decltype(syrk(alpha, a, Ret({size(a), size(a)}, get_allocator(a))))> {
	return syrk(alpha, a, Ret({size(a), size(a)}, get_allocator(a)));  }

template<class A2D>
[[nodiscard]]
auto syrk(A2D const& A)  // NOLINT(readability-identifier-length) BLAS naming
->decltype(syrk(1., A)) {
	return syrk(1., A); }

} // end namespace boost::multi::blas

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS syrk"
//#include<boost/test/unit_test.hpp>

//#include "../blas/gemm.hpp"

//#include "../../array.hpp"
//#include "../../utility.hpp"

//#include <boost/timer/timer.hpp>

//#include<complex>
//#include<cassert>
//#include<iostream>
//#include<numeric>
//#include<algorithm>

////#include<catch.hpp>

//using std::cout;
//using std::cerr;

//namespace multi = boost::multi;

//template<class M> decltype(auto) print(M const& C){
//	using boost::multi::size;
//	for(int i = 0; i != size(C); ++i){
//		for(int j = 0; j != size(C[i]); ++j)
//			std::cout << C[i][j] << ' ';
//		std::cout << std::endl;
//	}
//	return std::cout << std::endl;
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_real){
//	multi::array<double, 2> const a = {
//		{ 1., 3., 4.},
//		{ 9., 7., 1.}
//	};
//	{
//		multi::array<double, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		using blas::transposed;
//		syrk(filling::lower, 1., transposed(a), 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[2][1] == 19. ); 
//		BOOST_REQUIRE( c[1][2] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		using blas::transposed;
//		syrk(filling::upper, 1., transposed(a), 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[1][2] == 19. );
//		BOOST_REQUIRE( c[2][1] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		using blas::syrk;
//		syrk(filling::lower, 1., a, 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0] == 34. ); 
//		BOOST_REQUIRE( c[0][1] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		syrk(filling::upper, 1., a, 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, a⸆a, `c` in lower triangular
//		BOOST_REQUIRE( c[0][1] == 34. ); 
//		BOOST_REQUIRE( c[1][0] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		syrk(filling::upper, 1., a, 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, a⸆a, `c` in lower triangular
//		BOOST_REQUIRE( c[0][1] == 34. ); 
//		BOOST_REQUIRE( c[1][0] == 9999. );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_real_special_case){
//	multi::array<double, 2> const a = {
//		{ 1., 3., 4.},
//	};
//	{
//		multi::array<double, 2> c({1, 1}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		syrk(filling::lower, 1., a, 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		//BOOST_REQUIRE( c[1][0] == 34. ); 
//		//BOOST_REQUIRE( c[0][1] == 9999. );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_complex_real_case){
//	using complex = std::complex<double>;
//	multi::array<complex, 2> const a = {
//		{ 1., 3., 4.},
//		{ 9., 7., 1.}
//	};
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		using blas::transposed;
//		syrk(filling::lower, 1., transposed(a), 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[2][1] == 19. );
//		BOOST_REQUIRE( c[1][2] == 9999. );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_complex){
//	using complex = std::complex<double>;
//	constexpr auto const I = complex{0., 1.};
//	multi::array<complex, 2> const a = {
//		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
//		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
//	};
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		using blas::transposed;
//		syrk(filling::lower, 1., transposed(a), 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[2][1] == complex(-3., -34.) );
//		BOOST_REQUIRE( c[1][2] == 9999. );
//	}
//	{
//		multi::array<complex, 2> c({2, 2}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		syrk(filling::lower, 1., a, 0., c); // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0] == complex(18., -21.) );
//		BOOST_REQUIRE( c[0][1] == 9999. );
//	}
//	{
//		multi::array<complex, 2> c({2, 2}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		syrk(filling::upper, 1., a, 0., c); // c⸆=c=aa⸆=(aa⸆)⸆, `c` in upper triangular
//		BOOST_REQUIRE( c[0][1] == complex(18., -21.) ); 
//		BOOST_REQUIRE( c[1][0] == 9999. );
//	}
//}


//BOOST_AUTO_TEST_CASE(multi_blas_syrk_automatic_operation_complex){
//	using complex = std::complex<double>;
//	constexpr auto const I = complex{0., 1.};
//	multi::array<complex, 2> const a = {
//		{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
//		{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
//	};
//	{
//		multi::array<complex, 2> c({2, 2}, 9999.);
//		using multi::blas::filling;
//		syrk(filling::lower, 1., a, 0., c); // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0]==complex(18., -21.) );
//		BOOST_REQUIRE( c[0][1]==9999. );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		using blas::transposed;
//		syrk(filling::lower, 1., transposed(a), 0., c); // c⸆=c=a⸆a=(aa⸆)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[2][1]==complex(-3.,-34.) );
//		BOOST_REQUIRE( c[1][2]==9999. );
//	}
//	{
//		multi::array<complex, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		using blas::transposed;
//		syrk(filling::lower, 1., rotated(a), 0., c); // c⸆=c=a⸆a=(aa⸆)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[2][1]==complex(-3.,-34.) );
//		BOOST_REQUIRE( c[1][2]==9999. );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_automatic_operation_real){
//	multi::array<double, 2> const a = {
//		{ 1., 3., 4.},
//		{ 9., 7., 1.}
//	};
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		using multi::blas::filling;
//		syrk(filling::lower, 1., a, 0., c); // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0] == 34. );
//		BOOST_REQUIRE( c[0][1] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		using multi::blas::filling;
//		syrk(filling::upper, 1., a, 0., c); // c⸆=c=aa⸆=(aa⸆)⸆, `c` in upper triangular
//		BOOST_REQUIRE( c[0][1] == 34. );
//		BOOST_REQUIRE( c[1][0] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({3, 3}, 9999.);
//		using multi::blas::filling;
//		syrk(filling::lower, 1., rotated(a), 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[2][1] == 19. );
//		BOOST_REQUIRE( c[1][2] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::transposed;
//		using blas::filling;
//		syrk(filling::lower, 1., transposed(a), 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[2][1] == 19. );
//		BOOST_REQUIRE( c[1][2] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({3, 3}, 9999.);
//		namespace blas = multi::blas;
//		using blas::transposed;
//		using blas::filling;
//		syrk(filling::upper, 1., transposed(a), 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in upper triangular
//		BOOST_REQUIRE( c[1][2] == 19. );
//		BOOST_REQUIRE( c[2][1] == 9999. );
//	}
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		using multi::blas::filling;
//		using multi::blas::transposed;
//		syrk(filling::upper, 1., a, 0., transposed(c)); // c⸆=c=aa⸆=(aa⸆)⸆, `c` in upper triangular
//		BOOST_REQUIRE( c[0][1] == 9999. ); 
//		BOOST_REQUIRE( c[1][0] == 34. );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_automatic_implicit_zero){
//	multi::array<double, 2> const a = {
//		{ 1., 3., 4.},
//		{ 9., 7., 1.}
//	};
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		using multi::blas::filling;
//		syrk(filling::lower, 1., a, c); // c⸆=c=aa⸆=(aa⸆)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0] == 34. );
//		BOOST_REQUIRE( c[0][1] == 9999. );
//	}
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_automatic_symmetrization){
//	multi::array<double, 2> const a = {
//		{ 1., 3., 4.},
//		{ 9., 7., 1.}
//	};
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		using multi::blas::syrk;
//		using multi::blas::gemm;
//		using multi::blas::T;
//		syrk(1., a, c); // c⸆=c=aa⸆=(aa⸆)⸆
//		BOOST_REQUIRE( c[1][0] == 34. );
//		BOOST_REQUIRE( c[0][1] == 34. );
//		BOOST_REQUIRE( syrk(a) == gemm(a, T(a)) );
//	}
//	{
//		using multi::blas::syrk;
//		multi::array<double, 2> c = syrk(1., a); // c⸆=c=aa⸆=(aa⸆)⸆
//		BOOST_REQUIRE( c[1][0] == 34. );
//		BOOST_REQUIRE( c[0][1] == 34. );
//	}
//	{
//		using multi::blas::syrk;
//		multi::array<double, 2> c = syrk(a); // c⸆=c=aa⸆=(aa⸆)⸆
//		BOOST_REQUIRE( c[1][0] == 34. );
//		BOOST_REQUIRE( c[0][1] == 34. );
//	}
//	{
//		using multi::blas::transposed;
//		using multi::blas::syrk;
//		multi::array<double, 2> c = syrk(transposed(a)); // c⸆=c=a⸆a=(a⸆a)⸆
//		BOOST_REQUIRE( c[2][1] == 19. );
//		BOOST_REQUIRE( c[1][2] == 19. );
//	}
//}

//#if 0



//}







//}







//#if 0
//	{

//		{
//			multi::array<complex, 2> C({2, 2}, 9999.);
//			syrk(1., rotated(A), rotated(C)); // C^T=C=A*A^T=(A*A^T)^T
//			assert( C[1][0] == complex(18., -21.) );
//		}
//		{
//			multi::array<complex, 2> C({2, 2}, 9999.);
//			syrk(rotated(A), rotated(C)); // C^T=C=A*A^T=(A*A^T)^T
//			assert( C[1][0] == complex(18., -21.) );
//		}
//		{
//			complex C[2][2];
//			using multi::rotated;
//			syrk(rotated(A), rotated(C)); // C^T=C=A*A^T=(A*A^T)^T
//			assert( C[1][0] == complex(18., -21.) );
//		}
//		{
//			auto C = syrk(1., A); // C = C^T = A^T*A, C is a value type matrix (with C-ordering, information is everywhere)
//			assert( C[1][2]==complex(-3.,-34.) );
//		}
//		{
////			what(rotated(syrk(A)));
//			multi::array C = rotated(syrk(A)); // C = C^T = A^T*A, C is a value type matrix (with C-ordering, information is in upper triangular part)
//			print(C) <<"---\n";
//		}
//		
//	}
//#if 0
//	{
//		multi::array<complex, 2> const A = {
//			{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
//			{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
//		};
//		auto C = rotated(syrk(A)).decay(); // C = C^T = A^T*A, C is a value type matrix (with C-ordering, information is in upper triangular part)
//		print(C) <<"---\n";
////		print(C) <<"---\n";
//	}
//	return 0;
//	{
//		multi::array<complex, 2> const A = {
//			{ 1. + 3.*I, 3.- 2.*I, 4.+ 1.*I},
//			{ 9. + 1.*I, 7.- 8.*I, 1.- 3.*I}
//		};
//		auto C = syrk(rotated(A)); // C = C^T = A^T*A, C is a value type matrix (with C-ordering)
//		print(C) <<"---\n";
//	}
//#endif
//#endif
//}

//BOOST_AUTO_TEST_CASE(multi_blas_syrk_herk_fallback){
//	multi::array<double, 2> const a = {
//		{ 1., 3., 4.},
//		{ 9., 7., 1.}
//	};
//	{
//		multi::array<double, 2> c({2, 2}, 9999.);
//		namespace blas = multi::blas;
//		using blas::filling;
//		syrk(filling::lower, 1., a, 0., c); // c⸆=c=a⸆a=(a⸆a)⸆, `c` in lower triangular
//		BOOST_REQUIRE( c[1][0] == 34. ); 
//		BOOST_REQUIRE( c[0][1] == 9999. );
//	}
//}
//#endif

//#endif
#endif

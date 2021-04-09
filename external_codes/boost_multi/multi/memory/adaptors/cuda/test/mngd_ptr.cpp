
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../../adaptors/cuda/cublas/context.hpp"

#include "../managed/ptr.hpp"
#include "../../../../adaptors/cuda.hpp"
#include "../../../../adaptors/blas/gemm.hpp"
#include "../../../../adaptors/blas/trsm.hpp"

#include<random>

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;
namespace blas = multi::blas;

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(multi_cuda_mngd_ptr){
	using T = double;
	static_assert( sizeof(cuda::managed::ptr<T>) == sizeof(T*) );
	static_assert( std::is_convertible<cuda::managed::ptr<T>, T*>{} );
	auto f = [](double* dp){return bool{dp};};
	cuda::managed::ptr<T> p;
	f(p);
}

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_cuda_mngd_ptr_call_gemm){
	using complex = std::complex<double>; complex const I{0, 1};
	boost::multi::cuda::managed::array<complex, 2> m = {
		{ 1. + 2.*I, 3. - 3.*I, 1.-9.*I},
		{ 9. + 1.*I, 7. + 4.*I, 1.-8.*I},
	};
	boost::multi::cuda::managed::array<complex, 2> const b = {
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};
//	{
//		blas::context ctxt;
//		auto c =+ blas::gemm(&ctxt, 1., m, b);
//		static_assert( std::is_same<decltype(c), multi::cuda::managed::array<complex, 2>>{} );
//		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
//		BOOST_REQUIRE( b[1][2] == 2.+1.*I );
//	}
//	{
//		multi::cuda::managed::array<complex, 2> c({2, 4});
//		multi::cuda::cublas::context ctxt;
//		blas::gemm_n(ctxt, 1., begin(m), size(m), begin(b), 0., begin(c));
//		cudaDeviceSynchronize();
//		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
//		BOOST_REQUIRE( b[1][2] == 2.+1.*I );
//	}
	{
		multi::cuda::cublas::context ctxt;
		auto c =+ blas::gemm(&ctxt, 1., m, b);
		static_assert( std::is_same<decltype(c), multi::cuda::managed::array<complex, 2>>{} );
		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
		BOOST_REQUIRE( b[1][2] == 2.+1.*I );
	}
//	{
//		auto c =+ blas::gemm(1., m, b);
//		static_assert( std::is_same<decltype(c), multi::cuda::managed::array<complex, 2>>{} );
//		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
//		BOOST_REQUIRE( b[1][2] == 2.+1.*I );
//	}
//	{
//		auto c =+ blas::gemm(1., m, b);//blas::default_context_of(m.base()), 1., m, b);
//		static_assert( std::is_same<decltype(c), multi::cuda::managed::array<complex, 2>>{} );
//		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
//		BOOST_REQUIRE( b[1][2] == 2.+1.*I );
//	}
//	{
//		multi::cuda::managed::array<complex, 2> c({2, 4});
//		multi::cuda::cublas::context ctxt;
//		blas::gemm_n(ctxt, 1., begin(m), size(m), begin(b), 0., begin(c));
//		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
//		BOOST_REQUIRE( b[1][2] == 2.+1.*I );
//	}
//	{
//		auto c =+ blas::gemm(1., m, b);
//		static_assert( std::is_same<decltype(c), multi::cuda::managed::array<complex, 2>>{} );
//		BOOST_REQUIRE( c[1][2] == complex(112, 12) );
//		BOOST_REQUIRE( b[1][2] == 2.+1.*I );
//	}
}

//BOOST_AUTO_TEST_CASE(multi_blas_trsm_complex_nonsquare_default_diagonal_hermitized_gemm_check_managed_ptr, *utf::tolerance(0.00001)){
//	using complex = std::complex<double>; complex const I{0, 1};
//	multi::cuda::managed::array<complex, 2> const A = {
//		{ 1. + 4.*I,  3.,  4.- 10.*I},
//		{ 0.,  7.- 3.*I,  1.},
//		{ 0.,  0.,  8.- 2.*I}
//	};
//	namespace blas = multi::blas;
//	{
//		{
//			multi::cuda::managed::array<complex, 2> B = {
//				{1. + 1.*I, 5. + 3.*I},
//				{2. + 1.*I, 9. + 3.*I},
//				{3. + 1.*I, 1. - 1.*I},
//			};
//			blas::trsm(blas::side::left, blas::filling::lower, 1., blas::H(A), B); // S = A⁻¹†.B, S† = B†.A⁻¹
//			cudaDeviceSynchronize();
//			BOOST_TEST( real(B[2][1]) == 1.71608  );
//		}
//		{
//			multi::cuda::managed::array<complex, 2> B = {
//				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
//				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
//			};
//			auto const S =+ blas::trsm(blas::side::left, blas::filling::upper, 1., A, blas::H(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
//			cudaDeviceSynchronize();
//		//	BOOST_TEST( imag(S[2][1]) == +0.147059 );
//			BOOST_TEST( imag(B[1][2]) == -0.147059 );
//		}
//		{
//			multi::cuda::managed::array<complex, 2> B = {
//				{1. + 1.*I, 2. + 1.*I, 3. + 1.*I},
//				{5. + 3.*I, 9. + 3.*I, 1. - 1.*I}
//			};
//			auto const S =+ blas::trsm(blas::side::left, blas::filling::upper, 2., A, blas::H(B)); // S = A⁻¹B†, S†=B.A⁻¹†, S=(B.A⁻¹)†, B <- S†, B <- B.A⁻¹†
//			cudaDeviceSynchronize();
//		//	BOOST_TEST( imag(S[2][1]) == +0.147059*2. );
//			BOOST_TEST( imag(B[1][2]) == -0.147059*2. );
//		}
//	}
//}


#ifdef COMPILATION_INSTRUCTIONS
/usr/local/cuda-11.1/bin/nvcc -x cu -std=c++17  -use_fast_math -lpthread -D_REENTRANT -DBOOST_PP_VARIADICS  -Xcudafe "--diag_suppress=implicit_return_from_non_void_function" --extended-lambda --expt-relaxed-constexpr $0 -o $0x `pkg-config --cflags --libs cudart-11.0 cublas-11.0 blas` -lboost_unit_test_framework -DBOOST_LOG_DYN_LINK -lboost_log -lboost_thread -lboost_system -lboost_log_setup -lpthread -lboost_timer&&$0x&&rm $0x; exit
#endif
// © Alfredo A. Correa 2020-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi cuBLAS gemm"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

//#include"boost/log/trivial.hpp"
//#define MULTI_MARK_SCOPE(MsG) BOOST_LOG_TRIVIAL(trace)<<MsG

//#include "../../../../adaptors/cublas/context.hpp"

#include "../../../cuda/cublas.hpp"
#include "../../../../array.hpp"

#include "../../../../adaptors/cuda.hpp"
#include "../../../../adaptors/blas.hpp"

#include<random>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_gemm_complex_3x2_3x2){
	using complex = std::complex<double>; complex const I{0, 1};
	namespace blas = multi::blas;
	multi::array<complex, 2> const a = {
		{1. + 2.*I, 5. + 2.*I}, 
		{9. - 1.*I, 9. + 1.*I}, 
		{1. + 1.*I, 2. + 2.*I}
	};
	multi::array<complex, 2> const b = {
		{ 11. - 2.*I, 5. + 2.*I},
		{  7. - 3.*I, 2. + 1.*I},
		{  8. - 1.*I, 1. + 1.*I}
	};
	{
		{
			multi::array<complex, 2> c({2, 2});
			c = blas::gemm(1., blas::H(a), b); // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE( c[1][0] == 125.-84.*I );
		}
	}
	{
		multi::cuda::array<complex, 2> const a_gpu = a;
		multi::cuda::array<complex, 2> const b_gpu = b;
		{
			multi::cuda::array<complex, 2> c_gpu({2, 2});
			c_gpu = blas::gemm(1., blas::H(a_gpu), b_gpu); // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
		}
		{
			auto c_gpu =+ blas::gemm(1.0, blas::H(a_gpu), b_gpu);
			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
		}
	}
	{
		multi::cuda::managed::array<complex, 2> const a_gpu = a;
		multi::cuda::managed::array<complex, 2> const b_gpu = b;
		{
			multi::cuda::managed::array<complex, 2> c_gpu({2, 2});
			blas::gemm(1., blas::H(a_gpu), b_gpu, 0., c_gpu); // c=ab, c⸆=b⸆a⸆
			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
		}
		{
			auto c_gpu =+ blas::gemm(1.0, blas::H(a_gpu), b_gpu);
			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
		}
	}
}

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_gemm_complex_3x2_3x2_with_context){
//	using complex = std::complex<double>; complex const I{0, 1};
//	namespace blas = multi::blas;
//	multi::array<complex, 2> const a = {
//		{1. + 2.*I, 5. + 2.*I}, 
//		{9. - 1.*I, 9. + 1.*I}, 
//		{1. + 1.*I, 2. + 2.*I}
//	};
//	multi::array<complex, 2> const b = {
//		{ 11. - 2.*I, 5. + 2.*I},
//		{  7. - 3.*I, 2. + 1.*I},
//		{  8. - 1.*I, 1. + 1.*I}
//	};
//	{
//		{
//			multi::blas::context ctx;
//			multi::array<complex, 2> c({2, 2});
//			blas::gemm(ctx, 1., blas::H(a), b, 0., c); // c=ab, c⸆=b⸆a⸆
//			BOOST_REQUIRE( c[1][0] == 125.-84.*I );
//		}
//	}
//	{
//		multi::cublas::context ctx;
//		multi::cuda::array<complex, 2> const a_gpu = a;
//		multi::cuda::array<complex, 2> const b_gpu = b;
//		{
//			multi::cuda::array<complex, 2> c_gpu({2, 2});
//			blas::gemm(ctx, 1., blas::H(a_gpu), b_gpu, 0., c_gpu); // c=ab, c⸆=b⸆a⸆
//			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//		}
//		{
//			auto c_gpu =+ blas::gemm(&ctx, blas::H(a_gpu), b_gpu);
//			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//		}
//	}
//	{
//		multi::cublas::context ctx;
//		multi::cuda::managed::array<complex, 2> const a_gpu = a;
//		multi::cuda::managed::array<complex, 2> const b_gpu = b;
//		{
//			multi::cuda::managed::array<complex, 2> c_gpu({2, 2});
//			blas::gemm(ctx, 1., blas::H(a_gpu), b_gpu, 0., c_gpu); // c=ab, c⸆=b⸆a⸆
//			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//		}
//		{
//			auto c_gpu =+ blas::gemm(&ctx, blas::H(a_gpu), b_gpu);
//			BOOST_REQUIRE( c_gpu[1][0] == 125.-84.*I );
//		}
//	}
//}

#if 0
BOOST_AUTO_TEST_CASE(multi_adaptors_blas_cuda_gemm_context_timing){
	using complex = std::complex<double>;//complex const I{0, 1};
	
	multi::array<complex, 2> A({1000, 1000});
	multi::array<complex, 2> B(      {1000, 1000});
	multi::array<complex, 2> C({size(A), size(~B)});
	A[99][99] = B[11][22] = C[33][44] = 1.0;
	std::cerr<< "memory " << (A.num_elements()+ B.num_elements() + C.num_elements())*sizeof(complex)/1e6 <<" MB"<<std::endl;
	
	{
		auto rand = [d=std::uniform_real_distribution<>{0., 10.}, g=std::mt19937{}]() mutable{return complex{d(g), d(g)};};
		std::generate(A.elements().begin(), A.elements().end(), rand);
		std::generate(B.elements().begin(), B.elements().end(), rand);
	}
	namespace blas = multi::blas;
	{
		boost::timer::auto_cpu_timer t; // 2.398206s
		for(auto i = 0; i != 10; ++i){
			blas::context ctx;
			blas::gemm(ctx, 1, A, B, 0, C);
		}
	}
	using device_array = multi::cuda::array<complex, 2>;
	{
		device_array A_gpu = A, B_gpu = B, C_gpu({size(A), size(~B)});

		boost::timer::auto_cpu_timer t; // 0.707426s
		for(auto i = 0; i != 10; ++i){
			multi::cublas::context ctx;
			blas::gemm(ctx, 1, A_gpu, B_gpu, 0, C_gpu);
		}
	}
	{
		device_array A_gpu = A, B_gpu = B, C_gpu({size(A), size(~B)});

		boost::timer::auto_cpu_timer t; // 0.613534s
		multi::cublas::context ctx;
		for(auto i = 0; i != 10; ++i) blas::gemm(ctx, 1, A_gpu, B_gpu, 0, C_gpu);
	}
}
#endif


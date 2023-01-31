#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
nvcc    -D_TEST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY -x cu                                    $0 -o $0x          -lboost_unit_test_framework -lboost_timer&&$0x&&
clang++ -D_TEST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY -x cuda --cuda-gpu-arch=sm_61 -std=c++14 $0 -o $0x -lcudart -lboost_unit_test_framework -lboost_timer&&$0x&&
rm $0x; exit
#endif

#ifndef MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY_HPP
#define MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY_HPP

#include<cassert>
//#include<iostream>

#include "../../../adaptors/cuda.hpp"
//#include "../algorithms/for_each.hpp"

//#include "/home/correaa/prj/alf/boost/iterator/zipper.hpp"

#ifndef HD
#if defined(__CUDACC__)
#define HD __host__ __device__
#else
#define HD
#endif
#endif

namespace boost{
namespace multi{namespace cuda{

#if 0
template<typename From, typename To, typename = std::enable_if_t<std::is_trivially_assignable<To&, From>{}> >
array_iterator<To, 1, To*> copy(
	array_iterator<From, 1, memory::cuda::ptr<To>> f, 
	array_iterator<From, 1, memory::cuda::ptr<To>> l, 
	array_iterator<To, 1, To*> d
){
	assert(0);
	assert(f.stride() == l.stride()); static_assert(sizeof(From) == sizeof(To), "!");
	auto n = std::distance(f, l);
	if(f.stride()==1 and d.stride()==1){
		auto s = cudaMemcpy(d.data(), raw_pointer_cast(f.data()), n*sizeof(To), cudaMemcpyDeviceToHost); assert( s == cudaSuccess );
	}else{
		auto s = cudaMemcpy2D(d.data(), d.stride()*sizeof(To), raw_pointer_cast(f.data()), f.stride()*sizeof(To), sizeof(To), n, cudaMemcpyDeviceToHost);
		assert( s == cudaSuccess );
	}
	return d + n;
}

template<typename From, typename From2, typename To, typename To2, typename = std::enable_if_t<std::is_trivially_assignable<To&, From>{}> >
array_iterator<To, 1, To*> copy(
	array_iterator<From, 1, memory::cuda::ptr<From2>> f, 
	array_iterator<From, 1, memory::cuda::ptr<From2>> l, 
	array_iterator<To  , 1, memory::cuda::ptr<To2>  > d
){
	assert(0);
	assert(f.stride() == l.stride()); static_assert(sizeof(From) == sizeof(To), "!");
	auto n = std::distance(f, l);
	if(f.stride()==1 and d.stride()==1){
		auto s = cudaMemcpy(raw_pointer_cast(d.data()), raw_pointer_cast(f.data()), n*sizeof(To), cudaMemcpyDeviceToHost); assert( s == cudaSuccess );
	}else{
		auto s = cudaMemcpy2D(raw_pointer_cast(d.data()), d.stride()*sizeof(To), raw_pointer_cast(f.data()), f.stride()*sizeof(To), sizeof(To), n, cudaMemcpyDeviceToDevice);
		assert( s == cudaSuccess );
	}
	return d + n;
}
#endif

}}
}


#ifdef _TEST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA copy"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../adaptors/cuda.hpp"

#include <thrust/for_each.h>
#include <thrust/execution_policy.h>

#include <boost/timer/timer.hpp>

#if __cpp_lib_parallel_algorithm >= 201603
#include<execution>
#endif

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

template<class T> __device__ void WHAT(T&&) = delete;
template<class T> __device__ void WHAT(int) = delete;

template<class T> T&& what(T&&) = delete;

BOOST_AUTO_TEST_CASE(copy_by_iterator){
	auto const A_cpu = []{
		multi::array<double, 2> r({198, 23});
		std::generate(r.data_elements(), r.data_elements()+r.num_elements(), &std::rand);
		return r;
	}();
	multi::cuda::array<double, 2> A = A_cpu;

	multi::cuda::array<double, 2> B(extensions(A));
	B() = A();
//	BOOST_REQUIRE( A[13] == B[13] );
}

BOOST_AUTO_TEST_CASE(copy_by_pointer){
	auto const A_cpu = []{
		multi::array<double, 2> r({198, 23});
		std::generate(r.data_elements(), r.data_elements()+r.num_elements(), &std::rand);
		return r;
	}();
	multi::cuda::array<double, 2> A = A_cpu;

	multi::cuda::array<double, 2> B(extensions(A));
	B = A;
//	BOOST_REQUIRE( A[13] == B[13] );
}


BOOST_AUTO_TEST_CASE(cuda_copy){

	multi::cuda::array<double, 1> A(1<<27); CUDA_SLOW( A[10] = 99. );
	multi::cuda::array<double, 1> B(size(A));

	{
		boost::timer::auto_cpu_timer t{"thrust copy_n cuda::ptr %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		thrust::copy_n(thrust::device, A.data_elements(), A.num_elements(), B.data_elements());
	}
	{
		boost::timer::auto_cpu_timer t{"cuda copy_n cuda::ptr copy_n %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		copy_n(A.data_elements(), A.num_elements(), B.data_elements());
	}
	{
		boost::timer::auto_cpu_timer t{"cuda copy_n cuda::ptr copy_n %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		multi::adl::copy_n(A.data_elements(), A.num_elements(), B.data_elements());
	}
#if 0
	{
		boost::timer::auto_cpu_timer t{"cuda ptr copy_n %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		/*multi::cuda::*/copy_n(A.data_elements(), A.num_elements(), B.data_elements());
	}
	{
		boost::timer::auto_cpu_timer t{"indirect cuda ptr copy_n %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		B = A;
	}
	{
		boost::timer::auto_cpu_timer t{"indirect cuda ptr uninitialized_copy_n %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		multi::cuda::array<double, 1> C = A;
		BOOST_REQUIRE( CUDA_SLOW( C[10] == 99. ) );
	}
	{
		boost::timer::auto_cpu_timer t{"indirect cuda ptr uninitialized_copy_n %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		multi::cuda::array<double, 1> C = A;//();
		BOOST_REQUIRE( CUDA_SLOW( C[10] == 99. ) );
	}
	BOOST_REQUIRE( CUDA_SLOW( B[10] == 99. ) );
	CUDA_SLOW( B[10] = 10. );
	{
		boost::timer::auto_cpu_timer t{"thrust copy_n %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		thrust::copy_n(thrust::device, begin(A), size(A), begin(B));
	}
	BOOST_REQUIRE( CUDA_SLOW( B[10] == 99. ) );
#endif

/*	multi::cuda::for_each_n(
		boost::iterators::zip(begin(A), begin(B)), 
		size(A), 
		[]__device__(auto&& e){
			std::get<1>(e) = std::get<0>(e);
			printf( "**** %f %f\n", static_cast<double const&>(std::get<0>(e)), static_cast<double const&>(std::get<1>(e)) ); 
		}
	);*/

//	auto l = 
//	BOOST_REQUIRE( l == end(B) );
//	std::cout << B[8] << std::endl;
//	multi::cuda::array<double, 1> A(10, 99.);
//	BOOST_REQUIRE( CUDA_SLOW( A[5] == 99. ) );
//	int uno = 1.;
//	for_each(begin(A), end(A), [uno]__device__(auto&& e){e = uno;});
//	BOOST_REQUIRE( CUDA_SLOW( A[5] == 1. ) );
}

#if 0
BOOST_AUTO_TEST_CASE(cuda_for_each){
	multi::cuda::array<double, 1> A(10, 99.);
	BOOST_REQUIRE( CUDA_SLOW( A[5] == 99. ) );
	int uno = 1.;
	for_each(begin(A), end(A), [uno]__device__(auto&& e){e = uno;});
	BOOST_REQUIRE( CUDA_SLOW( A[5] == 1. ) );
}

BOOST_AUTO_TEST_CASE(cuda_timing){
	multi::cuda::managed::array<double, 1> A(1<<29); //std::cout << A.size()*8 << std::endl; 
	{
		boost::timer::auto_cpu_timer t{"cuda cold %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		multi::cuda::for_each(begin(A), end(A), []__device__(auto&& e){e = 11.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 11.) );
	{
		boost::timer::auto_cpu_timer t{"cuda %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		multi::cuda::for_each(begin(A), end(A), []__device__(auto&& e){e = 22.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 22.) );
	{
		boost::timer::auto_cpu_timer t{"thrust %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		thrust::for_each(thrust::device, begin(A), end(A), []__device__(auto&& e){e = 222.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 222.) );
	{
		std::for_each(begin(A), end(A), [](auto&& e){e = 55.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 55.) );
#if __cpp_lib_parallel_algorithm >= 201603
	{
		boost::timer::auto_cpu_timer t{"par %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		std::for_each(std::execution::par_unseq, begin(A), end(A), [](auto&& e){e = 33.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 33.) );
#endif
	{
		boost::timer::auto_cpu_timer t{"seq %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		std::for_each(begin(A), end(A), [](auto&& e){e = 55.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 55.) );
	{
		boost::timer::auto_cpu_timer t{"cuda cold %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		multi::cuda::for_each(begin(A), end(A), []__device__(auto&& e){e = 66.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 66.) );
	{
		boost::timer::auto_cpu_timer t{"cuda %ws wall, %us user + %ss system = %ts CPU (%p%)\n"};
		multi::cuda::for_each(begin(A), end(A), []__device__(auto&& e){e = 77.;});
	} BOOST_REQUIRE( CUDA_SLOW( A[size(A) - 10] == 77.) );
}
#endif

#endif
#endif


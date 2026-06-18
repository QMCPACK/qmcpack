#ifdef COMPILATION_INSTRUCTIONS//-*-indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*-
nvcc    -D_TEST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY -x cu                                    $0 -o $0x          -lboost_unit_test_framework -lboost_timer&&$0x&&
clang++ -D_TEST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY -x cuda --cuda-gpu-arch=sm_61 -std=c++14 $0 -o $0x -lcudart -lboost_unit_test_framework -lboost_timer&&$0x&&
rm $0x; exit
#endif

#ifndef BOOST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY_HPP
#define BOOST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY_HPP

#include<cassert>

#include "../../../adaptors/cuda.hpp"
//#include "../algorithms/for_each.hpp"

#ifndef BOOST_MULTI_HD
#if defined(__CUDACC__)
#define BOOST_MULTI_HD __host__ __device__
#else
#define BOOST_MULTI_HD
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

#endif  // BOOST_MULTI_ADAPTORS_CUDA_ALGORITHMS_COPY_HPP

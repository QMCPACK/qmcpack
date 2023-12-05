#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
$CXXX $CXXFLAGS $0 -o $0x -lcudart -lboost_unit_test_framework -lboost_timer&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef BOOST_MULTI_MEMORY_ADAPTORS_CUDA_CSTRING_HPP
#define BOOST_MULTI_MEMORY_ADAPTORS_CUDA_CSTRING_HPP

#include "../../adaptors/cuda/ptr.hpp"
#include "../../adaptors/cuda/managed/ptr.hpp"

#include<cuda_runtime.h> // cudaMemcpy/cudaMemset

#include<iostream>

namespace boost{
namespace multi{
namespace memory{
namespace cuda{
#if (__cpp_nontype_template_parameter_auto>=201606) or defined(__NVCC__)
template<auto CudaFunction, class... Args> // requires c++17
void call(Args&&... args){
	auto s = static_cast<Cuda::error>(CudaFunction(args...));
	if( s != Cuda::error::success ) throw std::system_error{make_error_code(s), "cannot call cuda function "};
}
#endif

template<class T, T CudaFunction>
auto call_static(std::string const& name){
	return [=](auto... args)->decltype(CublasFunction(args...), void()){
		std::cerr << "Calling function " << name << std::endl;
		Cuda::error s = CublasFunction(args...);
		if( s != Cuda::error::success ) throw std::system_error{make_error_code(s), "cannot call cuda function "};
	};
}

//template<class T, T CudaFunction>
//auto call_static(){return call_static(std::string{});}

#define CUDA(FunctionPostfix) ::boost::multi::memory::cuda::call_static<decltype(&cuda##FunctionPostfix), cuda##FunctionPostfix>(#FunctionPostfix)

namespace memcpy_{
//https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__TYPES.html#group__CUDART__TYPES_1g18fa99055ee694244a270e4d5101e95b
	enum class kind : std::underlying_type_t<enum cudaMemcpyKind>{
		host_to_host=cudaMemcpyHostToHost, host_to_device=cudaMemcpyHostToDevice,
		device_to_host=cudaMemcpyDeviceToHost, device_to_device=cudaMemcpyDeviceToDevice,
		inferred = cudaMemcpyDefault, default_ = cudaMemcpyDefault
	};
	template<class T1, class T2> constexpr kind type(T1*    , T2*    ){return kind::host_to_host    ;}
	template<class T1, class T2> constexpr kind type(ptr<T1>, T2*    ){return kind::host_to_device  ;}
	template<class T1, class T2> constexpr kind type(T1*    , ptr<T2>){return kind::device_to_host  ;}
	template<class T1, class T2> constexpr kind type(ptr<T1>, ptr<T2>){return kind::device_to_device;}
	template<class T1, class P2> constexpr kind type(managed::ptr<T1>, P2){return kind::inferred;}
	template<class P1, class T2> constexpr kind type(P1, managed::ptr<T2>){return kind::inferred;}
	template<class T1, class T2> constexpr kind type(managed::ptr<T1>, managed::ptr<T2>){return kind::inferred;}
	[[deprecated]] constexpr kind type(...)        {return kind::inferred;        }
}

template<typename Dest, typename Src, typename = decltype(memcpy_::type(Dest{}, Src{}))>
Dest memcpy(Dest dest, Src src, std::size_t byte_count){
	cuda::call<cudaMemcpy>(static_cast<void*>(dest), static_cast<void const*>(src), byte_count, static_cast<cudaMemcpyKind>(memcpy_::type(dest, src)));
	return dest;
}

inline ptr<void> memset(ptr<void> dest, int ch, std::size_t byte_count){
	cuda::call<cudaMemset>(static_cast<void*>(dest), ch, byte_count);
	return dest;
}

template<class VoidPDst = void*, class VoidPCSrc = void const*>
auto memcpy2D(VoidPDst dst, std::size_t dpitch, VoidPCSrc src, std::size_t spitch, std::size_t width, std::size_t height)
->decltype(cuda::call<cudaMemcpy2D>(static_cast<void*>(dst), dpitch, static_cast<void const*>(src), spitch, width, height, static_cast<cudaMemcpyKind>(memcpy_::type(dst, src)))){
	return cuda::call<cudaMemcpy2D>(static_cast<void*>(dst), dpitch, static_cast<void const*>(src), spitch, width, height, static_cast<cudaMemcpyKind>(memcpy_::type(dst, src)));}

}}}}

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUDA cstring"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include "../../adaptors/cuda/allocator.hpp"

//#include<boost/timer/timer.hpp>

//#include<numeric>

//namespace multi = boost::multi;
//namespace cuda = multi::memory::cuda;

//BOOST_AUTO_TEST_CASE(multi_memory_cuda_cstring){

//	std::size_t const n = 2e9/sizeof(double);
//	cuda::ptr<double> p = cuda::allocator<double>{}.allocate(n);
//	{
//		boost::timer::auto_cpu_timer t;
//		memset(p, 0, n*sizeof(double));
//	}
//	BOOST_REQUIRE( p[n/2]==0 );
//	CUDA_SLOW ( 
//		p[n/2] = 99.;
//	)
//	cuda::ptr<double> q = cuda::allocator<double>{}.allocate(n);
//	{
//		boost::timer::auto_cpu_timer t;
//		memcpy(q, p, n*sizeof(double));
//	}
//	BOOST_REQUIRE( p[n/2] == 99. );
//	BOOST_REQUIRE( q[n/2] == 99. );

//	double a = 5.;
//	BOOST_REQUIRE(a == 5.);

//}
//#endif
#endif

#ifdef COMPILATION_INSTRUCTIONS
(echo '#include "'$0'"'>$0.cpp)&&nvcc -x cu `#-Wall -Wextra` -D_TEST_MULTI_MEMORY_ADAPTOR_CUDA_CLIB $0.cpp -o $0x -lcudart&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_MEMORY_ADAPTOR_CUDA_CLIB_HPP
#define MULTI_MEMORY_ADAPTOR_CUDA_CLIB_HPP

#include<cuda_runtime.h> // cudaMalloc

#include    "../../adaptors/cuda/error.hpp"
#include "../../../config/NODISCARD.hpp"

namespace Cuda{
	using namespace std::string_literals;

	using size_t = ::size_t;
	error Malloc(void** p, size_t bytes){return static_cast<error>(cudaMalloc(p, bytes));}
	NODISCARD("because it will produce a memory leak")
	void* malloc(size_t bytes){
		void* ret;
		switch(auto e = Malloc(&ret, bytes)){
			case success           : return ret;
			case memory_allocation : return nullptr;
			default                : 
				throw std::system_error{e, "cannot allocate "+std::to_string(bytes)+" bytes in '"+__PRETTY_FUNCTION__+"'"};
		}
	}
	error Free(void* p){return static_cast<error>(cudaFree(p));}
	void free(void* p){
		auto e = Free(p);
		// probably will terminate if called from noexcept functon
		if(Cuda::success!=e) throw std::system_error{e, "cannot "s +__PRETTY_FUNCTION__}; 
	}
	namespace pointer{
		using attributes_t = cudaPointerAttributes;
		error GetAttributes(attributes_t* ret, void* p){return static_cast<error>(cudaPointerGetAttributes(ret, p));}
		attributes_t attributes(void* p){
			attributes_t ret;
			auto e = GetAttributes(&ret, p);
			if(e!=success) throw std::system_error{e, "cannot "s+__PRETTY_FUNCTION__};
			return ret;
		}
		bool is_device(void* p){return attributes(p).devicePointer or p==nullptr;}
	}
}

#ifdef _TEST_MULTI_MEMORY_ADAPTOR_CUDA_CLIB

#include "../cuda/ptr.hpp"
#include "../cuda/cstring.hpp"

#include<iostream>

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

using std::cout;

int main(){
	{
		void* p = Cuda::malloc(100);
		Cuda::free(p);
	}
	{
		char* p = (char*)Cuda::malloc(1ul<<50);
		assert(!p);
		Cuda::free(p);
	}
}
#endif
#endif


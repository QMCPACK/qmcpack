//#ifdef COMPILATION_INSTRUCTIONS
//(echo '#include "'$0'"'>$0.cpp)&&c++ -std=c++11 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_MEMORY_ADAPTOR_CUDA_MANAGED_MALLOC $0.cpp -lcudart -o $0x &&$0x&& rm $0x $0.cpp; exit
//#endif
#ifndef MULTI_MEMORY_ADAPTOR_CUDA_MANAGED_CLIB_HPP
#define MULTI_MEMORY_ADAPTOR_CUDA_MANAGED_CLIB_HPP

#include<cuda_runtime.h>  // cudaMallocManaged

#include "../../../adaptors/cuda/clib.hpp"  // Cuda::free
#include "../../../adaptors/cuda/error.hpp"

namespace Cuda {
	namespace Managed {
		inline error Malloc(void** p, size_t bytes) {return static_cast<error>(cudaMallocManaged(p, bytes/*, cudaMemAttachGlobal*/));}
		inline void* malloc(size_t bytes) {
			void* ret;
			switch(auto e = Malloc(&ret, bytes)) {
				case success           : return ret;
				case memory_allocation : return nullptr;
				default                : 
					throw std::system_error{e, "cannot allocate "+std::to_string(bytes)+" bytes in '"+__PRETTY_FUNCTION__+"'"};
			}
		}
		inline void free(void* p) {return Cuda::free(p);}
	}
}

//#ifdef _TEST_MULTI_MEMORY_ADAPTOR_CUDA_MANAGED_MALLOC

//#include "../../cuda/managed/ptr.hpp"

//#include<iostream>

//namespace multi = boost::multi;
//namespace cuda = multi::memory::cuda;

//using std::cout;

//int main(){
//	void* p = Cuda::Managed::malloc(100);
//	Cuda::Managed::free(p);
//}
//#endif
#endif

//#ifdef COMPILATION_INSTRUCTIONS
//(echo '#include"'$0'" '>$0.cpp)&& `#nvcc -ccbin=cuda-`c++ -D_TEST_MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_MALLOC $0.cpp -o $0x -lcudart &&$0x&&rm $0x; exit
//#endif

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_MALLOC_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_MALLOC_HPP

#include "../../../adaptors/cuda/managed/clib.hpp"
#include "../../../adaptors/cuda/managed/ptr.hpp"

namespace boost {namespace multi {
namespace memory {

namespace cuda {

namespace managed {
	[[nodiscard]]
	inline managed::ptr<void> malloc(size_t bytes) {
	    MULTI_MARK_SCOPE("cuda::managed::malloc");
	    return managed::ptr<void>{Cuda::Managed::malloc(bytes)};
	}
	inline void free(managed::ptr<void> p) {
	    MULTI_MARK_SCOPE("cuda::managed::free");
	    Cuda::Managed::free(static_cast<void*>(p));
	}
}

}

}
}}

#endif

//#ifdef COMPILATION_INSTRUCTIONS
//(echo '#include"'$0'" '>$0.cpp)&& `#nvcc -ccbin=cuda-`c++ -D_TEST_MULTI_MEMORY_ADAPTORS_CUDA_CACHED_MALLOC $0.cpp -o $0x -lcudart &&$0x&&rm $0x; exit
//#endif

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_CACHED_MALLOC_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_CACHED_MALLOC_HPP

#include "../../../adaptors/cuda/cached/clib.hpp"
#include "../../../adaptors/cuda/cached/ptr.hpp"

namespace boost {namespace multi {
namespace memory {

namespace cuda {

namespace cached {
	[[nodiscard]]
	inline cached::ptr<void> malloc(size_t bytes) {
	    MULTI_MARK_SCOPE("cuda::cached::malloc");
	    return cached::ptr<void>{Cuda::Cached::malloc(bytes)};
	}
	inline void free(cached::ptr<void> p) {
	    MULTI_MARK_SCOPE("cuda::cached::free");
	    Cuda::Cached::free(static_cast<void*>(p));
	}
}

}

}
}}

#endif

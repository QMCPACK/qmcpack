#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&nvcc -D_TEST_MULTI_MEMORY_ADAPTORS_CUDA_MALLOC $0.cpp -o $0x &&$0x&&rm $0x; exit
#endif

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_MALLOC
#define MULTI_MEMORY_ADAPTORS_CUDA_MALLOC

#include "../../adaptors/cuda/clib.hpp"
#include "../../adaptors/cuda/ptr.hpp"

namespace boost{namespace multi{
namespace memory{

namespace cuda{
	using size_t = Cuda::size_t;
#if __cplusplus >= 201703L
#if __has_cpp_attribute(nodiscard) >= 201603L
	[[nodiscard]]
#endif
#endif
	inline auto malloc(size_t bytes) -> ptr<void>{return ptr<void>{Cuda::malloc(bytes)};}
	inline void free(ptr<void> p){Cuda::free(p);}
}

}
}}

#ifdef _TEST_MULTI_MEMORY_ADAPTORS_CUDA_MALLOC

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

int main(){
	using cuda::ptr;
	ptr<double> p = static_cast<ptr<double>>(cuda::malloc(100*sizeof(double)));
	p[10] = 99.;
	cuda::free(p);
}

#endif
#endif


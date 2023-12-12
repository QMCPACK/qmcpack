// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_MALLOC_
#define MULTI_MEMORY_ADAPTORS_CUDA_MALLOC_

#include "../../adaptors/cuda/clib.hpp"
#include "../../adaptors/cuda/ptr.hpp"

namespace boost {namespace multi {
namespace memory {

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

}  // namespace memory
}  // namespace multi
}  // namespace	boost

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
#endif  // MULTI_MEMORY_ADAPTORS_CUDA_MALLOC_


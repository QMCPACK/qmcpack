// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_MALLOC_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_MANAGED_MALLOC_HPP

#include "../../../adaptors/cuda/managed/clib.hpp"
#include "../../../adaptors/cuda/managed/ptr.hpp"

namespace boost::multi {
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
}  // end namespace managed
}  // end namespace cuda
}  // end namespace memory
}  // end namespace boost::multi

#endif

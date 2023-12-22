#ifndef MULTI_ADAPTORS_CUDA_CUBLAS_CALL_HPP
#define MULTI_ADAPTORS_CUDA_CUBLAS_CALL_HPP

#include "../cublas/error.hpp"

#if defined(__NVCC__)
#include<cuda_runtime.h> // cudaDeviceSynchronize
#else
#include<hip/hip_runtime.h> // cudaDeviceSynchronize
#endif

#if defined(__NVCC__)
#define hicup(name) cuda##name
#define HICUP(name) CU##name
#else
#define hicup(name) hip##name
#define HICUP(name) HIP##name
#endif

namespace boost{
namespace multi::cuda::cublas{

template<auto Function, class... Args> // needs C++17
void call(Args... args){
	auto e = static_cast<enum cublas::error>(Function(args...));
	if(e != cublas::error::success) { throw std::system_error{e, "cannot call function "+ std::string{__PRETTY_FUNCTION__}}; }
}

#define CUBLAS_(F) call<F>

}
}

#undef hicup
#undef HICUP
#endif


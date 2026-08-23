// Copyright 2020-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_CUDA_CUBLAS_CALL_HPP
#define BOOST_MULTI_ADAPTORS_CUDA_CUBLAS_CALL_HPP
#pragma once

#include "boost/multi/adaptors/cuda/cublas/error.hpp"

#if !defined(MULTI_USE_HIP)
#include<cuda_runtime.h>  // cudaDeviceSynchronize
#else
#include<hip/hip_runtime.h>  // cudaDeviceSynchronize
#endif

#if !defined(MULTI_USE_HIP)
#define hicup(name) cuda##name
#define HICUP(name) CU##name
#else
#define hicup(name) hip##name
#define HICUP(name) HIP##name
#endif

namespace boost::multi::cuda::cublas{

template<auto Function, class... Args> // needs C++17
void call(Args... args){
	auto e = static_cast<enum cublas::error>(Function(args...));
	if(e != cublas::error::success) { throw std::system_error{e, "cannot call function "+ std::string{__PRETTY_FUNCTION__}}; }
}

#define CUBLAS_(F) call<F>

}

#undef hicup
#undef HICUP
#endif  // BOOST_MULTI_ADAPTORS_CUDA_CUBLAS_CALL_HPP

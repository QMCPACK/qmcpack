//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CUDA_TYPE_MAPPING_HPP
#define QMCPLUSPLUS_CUDA_TYPE_MAPPING_HPP

#include <type_traits>
#include "type_traits/type_mapping.hpp"
#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuComplex.h>
#else
#include <hip/hip_complex.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif

namespace qmcplusplus
{

// This saves us writing specific overloads with reinterpret casts for different std::complex to cuComplex types.
template<typename T>
using CUDATypeMap =
    typename std::disjunction<OnTypesEqual<T, float, float>,
                              OnTypesEqual<T, double, double>,
                              OnTypesEqual<T, float*, float*>,
                              OnTypesEqual<T, double*, double*>,
                              OnTypesEqual<T, float**, float**>,
                              OnTypesEqual<T, double**, double**>,
                              OnTypesEqual<T, std::complex<double>, cuDoubleComplex>,
                              OnTypesEqual<T, std::complex<float>, cuComplex>,
                              OnTypesEqual<T, std::complex<double>*, cuDoubleComplex*>,
                              OnTypesEqual<T, std::complex<float>**, cuComplex**>,
                              OnTypesEqual<T, std::complex<double>**, cuDoubleComplex**>,
                              OnTypesEqual<T, std::complex<float>*, cuComplex*>,
                              OnTypesEqual<T, const std::complex<double>*, const cuDoubleComplex*>,
                              OnTypesEqual<T, const std::complex<float>*, const cuComplex*>,
                              OnTypesEqual<T, const std::complex<float>**, const cuComplex**>,
                              OnTypesEqual<T, const std::complex<double>**, const cuDoubleComplex**>,
                              OnTypesEqual<T, const std::complex<float>* const*, const cuComplex* const*>,
                              OnTypesEqual<T, const std::complex<double>* const*, const cuDoubleComplex* const*>,
                              default_type<void>>::type;

template<typename T>
CUDATypeMap<T> castCUDAType(T var)
{
  return reinterpret_cast<CUDATypeMap<T>>(var);
}

} // namespace qmcplusplus

#endif // QMCPLUSPLUS_CUDA_TYPE_MAPPING_HPP

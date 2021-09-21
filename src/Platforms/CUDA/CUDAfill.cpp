//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "CUDAfill.hpp"
#include <stdexcept>
#include "CUDAruntime.hpp"

namespace qmcplusplus
{
template<typename T>
void CUDAfill_n(T* ptr, size_t n, const T& value)
{
  if (value != T())
    throw std::runtime_error("CUDAfill_n doesn't support fill non T() values!");
  // setting 0 value on each byte should be 0 for int, float and double.
  cudaErrorCheck(cudaMemset(ptr, 0, n * sizeof(T)), "Memset failed in CUDAfill_n!");
}

template void CUDAfill_n<int>(int* ptr, size_t n, const int& value);
template void CUDAfill_n<size_t>(size_t* ptr, size_t n, const size_t& value);

template void CUDAfill_n<float>(float* ptr, size_t n, const float& value);
template void CUDAfill_n<double>(double* ptr, size_t n, const double& value);

template void CUDAfill_n<std::complex<float>>(std::complex<float>* ptr, size_t n, const std::complex<float>& value);
template void CUDAfill_n<std::complex<double>>(std::complex<double>* ptr, size_t n, const std::complex<double>& value);
} // namespace qmcplusplus

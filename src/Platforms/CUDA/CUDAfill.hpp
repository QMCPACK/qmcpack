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


#ifndef QMCPLUSPLUS_CUDAFILL_H
#define QMCPLUSPLUS_CUDAFILL_H

#include <complex>
#include <stdexcept>
#include <cstddef>

namespace qmcplusplus
{
template<typename T>
void CUDAfill_n(T* ptr, size_t n, const T& value);

extern template void CUDAfill_n<int>(int* ptr, size_t n, const int& value);
extern template void CUDAfill_n<size_t>(size_t* ptr, size_t n, const size_t& value);

extern template void CUDAfill_n<float>(float* ptr, size_t n, const float& value);
extern template void CUDAfill_n<double>(double* ptr, size_t n, const double& value);

extern template void CUDAfill_n<std::complex<float>>(std::complex<float>* ptr,
                                                     size_t n,
                                                     const std::complex<float>& value);
extern template void CUDAfill_n<std::complex<double>>(std::complex<double>* ptr,
                                                      size_t n,
                                                      const std::complex<double>& value);

} // namespace qmcplusplus
#endif

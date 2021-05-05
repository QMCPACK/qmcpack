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

#ifndef QMCPLUSPLUS_TEST_KERNELS_HPP
#define QMCPLUSPLUS_TEST_KERNELS_HPP

#include <cuda_runtime_api.h>

namespace qmcplusplus
{
namespace testing
{

template<typename T>
cudaError_t checkValueCUDA(cudaStream_t hstream, T* device_value_ptr, T value, bool& result);

/** just an arbitrary struct for testing */
struct DualStruct
{
  int index;
  double value;
};
cudaError_t checkDualStruct(cudaStream_t hstream, DualStruct* device_struct_ptr, DualStruct dual_struct, bool& result);
  
extern template cudaError_t checkValueCUDA(cudaStream_t hstream, double* device_value_ptr, double value, bool& result);
extern template cudaError_t checkValueCUDA(cudaStream_t hstream, float* device_value_ptr, float value, bool& result);

} // namespace testing
} // namespace qmcplusplus

#endif

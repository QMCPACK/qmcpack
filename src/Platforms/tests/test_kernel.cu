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

#include "test_kernel.hpp"
#include "CUDA/CUDAallocator.hpp"

namespace qmcplusplus
{
namespace testing
{
template<typename T>
__global__ void checkValue_kernel(T* device_value_ptr, const T value, bool* result)
{
  if (*device_value_ptr == value)
    *result = true;
  else
    *result = false;
}

template<typename T>
cudaError_t checkValueCUDA(cudaStream_t hstream, T* device_value_ptr, T value, bool& result)
{
  CUDAAllocator<bool> bool_allocator;
  CUDAHostAllocator<bool> host_bool_allocator;
  bool* device_result = bool_allocator.allocate(1);
  bool* host_result   = host_bool_allocator.allocate(1);
  dim3 dim_block(1);
  dim3 dim_grid(1);
  checkValue_kernel<<<dim_grid, dim_block, 0, hstream>>>(device_value_ptr, value, device_result);
  cudaCheck(cudaStreamSynchronize(hstream));
  cudaError_t kernel_error = cudaPeekAtLastError();
  cudaCheck(cudaMemcpyAsync(host_result, device_result, sizeof(bool), cudaMemcpyDeviceToHost, hstream));
  cudaCheck(cudaStreamSynchronize(hstream));
  result = *host_result;
  bool_allocator.deallocate(device_result, 1);
  host_bool_allocator.deallocate(host_result, 1);
  return kernel_error;
}

template cudaError_t checkValueCUDA(cudaStream_t hstream, double* device_value_ptr, double value, bool& result);

} // namespace testing
} // namespace qmcplusplus

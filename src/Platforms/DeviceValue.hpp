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

#ifndef QMCPLUSPLUS_DEVICE_VALUE_HPP
#define QMCPLUSPLUS_DEVICE_VALUE_HPP
#include <memory>
#include <cuda_runtime_api.h>

namespace qmcplusplus
{
/** For single values that need to be device side.
 *  I think this is useful for values such as alpha's and beta's and constants that
 *  need to be available outside of a particular kernel.
 *  Especially for non batched cuBLAS calls this makes for clearer source.
 */
template<typename T, class DUAL_ALLOC>
class DeviceValue
{
private:
  DUAL_ALLOC allocator_;
  T* value_ptr_;
public:
  DeviceValue() = delete;
  DeviceValue(const DeviceValue&) = delete;
  DeviceValue(const DeviceValue&&) = delete;
  /** Constructor that makes sense to use.
   *  since we may not have a hstream at this point these are synchronous copies.
   *  So don't construct or copy objects containing DeviceValues during the course of a run!
   */
  DeviceValue(const T& t_val) : value_ptr_(nullptr)
  {
    static_assert(qmc_allocator_traits<DUAL_ALLOC>::is_dual_space, "DeviceValue intended to be used with dual allocator");
    value_ptr_  = std::allocator_traits<DUAL_ALLOC>::allocate(allocator_, 1);
    *value_ptr_ = t_val;
    cudaErrorCheck(cudaMemcpy(allocator_.get_device_ptr(), value_ptr_, sizeof(T), cudaMemcpyHostToDevice),
                   "failed to copy to device initial value of DeviceValue");
  }
  // Not sure this is correct
  T operator=(const T& t_val) { *value_ptr_ = t_val; }
  ~DeviceValue() { std::allocator_traits<DUAL_ALLOC>::deallocate(allocator_, value_ptr_, 1); }
  T get_value(){return *value_ptr_;}
  T* getPtr(){return value_ptr_;}
  T* getDevicePtr() { return allocator_.get_device_ptr(); }
  const T* getDevicePtr() const { return allocator_.get_device_ptr(); }
};
  
} // namespace qmcplusplus

#endif

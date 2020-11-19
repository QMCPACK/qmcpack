//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef HOST_BUFFER_MANAGER_HPP
#define HOST_BUFFER_MANAGER_HPP

#include <memory>
#include <cstddef>

#include "AFQMC/Memory/buffer_allocators.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// Class that manages the memory resource that generates allcators for host memory.
// Follows a monostate-type pattern. All variables are static and refer to a global instance
// of the resource.
class HostBufferManager
{
public:
  using generator_t = BufferAllocatorGenerator<host_memory_resource, host_constructor<char>>;

  template<class T>
  using allocator_t = typename generator_t::template allocator<T>;

  HostBufferManager(size_t size)
  {
    if (not generator)
      generator = std::make_unique<generator_t>(host_memory_resource{}, size, host_constructor<char>{});
  }

  HostBufferManager() { require(true); }

  ~HostBufferManager() {}

  void release()
  {
    if (generator)
    {
      generator.reset(nullptr);
      generator                    = nullptr;
      initialized_by_derived_class = false;
    }
  }

  generator_t& get_generator()
  {
    require(true);
    return *generator;
  }

protected:
  // protected constructor for use by derived class
  // this is done to allow for double initialization
  HostBufferManager(size_t size, bool)
  {
    require(true);
    if (initialized_by_derived_class)
      throw std::runtime_error("Error: Incorrect global state in protected constructor.\n");
    initialized_by_derived_class = true;
  }

  // static pointers to global objects
  //static generator_t* generator;
  static std::unique_ptr<generator_t> generator;

  static bool initialized_by_derived_class;

  void require(bool ini)
  {
    if (ini && not generator)
    {
      throw std::runtime_error("Error: Incorrect global state in require (found uninitialized).\n");
    }
    else if (not ini && generator)
    {
      throw std::runtime_error("Error: Incorrect global state in require (found initialized).\n");
    }
  }
};

} // namespace afqmc
} // namespace qmcplusplus
#endif

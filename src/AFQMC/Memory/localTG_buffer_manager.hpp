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

#ifndef LOCALTG_BUFFER_MANAGER_HPP
#define LOCALTG_BUFFER_MANAGER_HPP

#include <memory>
#include <cstddef>

#include "AFQMC/Memory/buffer_allocators.hpp"
#include "AFQMC/Memory/device_buffer_manager.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// Class that manages the memory resource that generates allocators for shared memory.
// Reserved for buffers with the localTG communicator.
// Follows a monostate-type pattern. All variables are static and refer to a global instance
// of the resource.
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
//  using DeviceBufferManager = LocalTGBufferManager;
class LocalTGBufferManager : public DeviceBufferManager
{
public:
  using DeviceBufferManager::allocator_t;
  using DeviceBufferManager::generator_t;

  LocalTGBufferManager(mpi3::shared_communicator& local, size_t size) : DeviceBufferManager(size, true) {}
  LocalTGBufferManager() : DeviceBufferManager() {}
  ~LocalTGBufferManager(){};

  void release() { DeviceBufferManager::release(); }

  generator_t& get_generator()
  {
    if (not initialized_by_derived_class)
      throw std::runtime_error("Error: Incorrect derived global state in get_generator. \n");
    return DeviceBufferManager::get_generator();
  }

private:
  using DeviceBufferManager::initialized_by_derived_class;
};
#else
class LocalTGBufferManager
{
public:
  using generator_t = BufferAllocatorGenerator<shm_memory_resource, shm_constructor<char>>;

  template<class T>
  using allocator_t = typename generator_t::template allocator<T>;

  LocalTGBufferManager(mpi3::shared_communicator& local, size_t size)
  {
    require(false);
    generator = std::make_unique<generator_t>(::shm::memory_resource_shm_ptr_with_raw_ptr_dispatch{&local}, size,
                                              shm_constructor<char>{local});
  }

  LocalTGBufferManager() { require(true); }

  ~LocalTGBufferManager() {}

  void release()
  {
    if (generator)
    {
      generator.reset(nullptr);
      generator = nullptr;
    }
  }

  generator_t& get_generator()
  {
    require(true);
    return *generator;
  }

protected:
  // static pointers to global objects
  //static generator_t* generator;
  static std::unique_ptr<generator_t> generator;

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
#endif

} // namespace afqmc
} // namespace qmcplusplus
#endif

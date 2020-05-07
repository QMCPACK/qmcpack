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

#include <memory>
#include "AFQMC/Memory/buffer_allocators.h"
#include "mpi3/shared_communicator.hpp"

namespace qmcplusplus {
namespace afqmc {

  auto host_buffer_generator = std::make_shared<host_allocator_generator_type> ( 
            host_memory_resource{},
            std::size_t(10*1024*1024),
            device_constructor<char>{}
        );
//#if defined(ENABLE_ACCELERATOR)
#if defined(ENABLE_CUDA)
  std::shared_ptr<device_allocator_generator_type> device_buffer_generator{nullptr};
#else  
// on cpu build, device and host are the same memory space
  std::shared_ptr<device_allocator_generator_type> device_buffer_generator=host_buffer_generator;
#endif
// on device, device and shm are the same memory space
// arch::INIT will set localTG_buffer_generator to device_buffer_generator
// on CPU build, localTG_buffer_generator will be built in GlobalTG 
  std::shared_ptr<localTG_allocator_generator_type> localTG_buffer_generator{nullptr};


  void update_buffer_generators() {
    if(host_buffer_generator == nullptr ) {
      app_error()<<" Uninitialized host buffer generator in update_buffer_generator."
                 <<std::endl;
      APP_ABORT("");
    }    
    host_buffer_generator->update();    
    if(device_buffer_generator!=nullptr) device_buffer_generator->update();    
    if(localTG_buffer_generator!=nullptr) localTG_buffer_generator->update();
  }  

  void make_localTG_buffer_generator(mpi3::shared_communicator& local, 
                                     std::size_t sz) {
#if defined(ENABLE_CUDA)
    // do nothing, since localTG_buffer_generator is in device and setup in arch::INIT
#else
    if(localTG_buffer_generator == nullptr)
      localTG_buffer_generator = std::make_shared<localTG_allocator_generator_type> (
                      shm::memory_resource_shm_ptr_with_raw_ptr_dispatch{local},
                      sz, device_constructor<char>{} );
#endif
  }

  void destroy_shm_buffer_generators() {
#if defined(ENABLE_CUDA)
    // do nothing 
#else
    if(localTG_buffer_generator != nullptr)
      localTG_buffer_generator.reset(); 
    localTG_buffer_generator = nullptr;
#endif
  }

}
}

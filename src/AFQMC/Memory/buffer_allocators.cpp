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

namespace qmcplusplus {
namespace afqmc {

  auto host_buffer_generator = std::make_shared<host_allocator_generator_type> ( 
            host_memory_resource{},
             std::size_t(10*1024*1024),
             std::allocator<char>{}
        );
//#if defined(ENABLE_ACCELERATOR)
#if defined(ENABLE_CUDA)
  std::shared_ptr<device_allocator_generator_type> device_buffer_generator{nullptr};
#else  
  std::shared_ptr<device_allocator_generator_type> device_buffer_generator=host_buffer_generator;
#endif

  void update_buffer_generators() {
    if(device_buffer_generator == nullptr or host_buffer_generator == nullptr ) {
      app_error()<<" Uninitialized buffer generators in update_buffer_generator."
                 <<std::endl;
      APP_ABORT("");
    }    
    host_buffer_generator->update();    
    if(device_buffer_generator!=nullptr) device_buffer_generator->update();    
  }  

/*
  void initialize_buffer_generators(std::size_t size_in_MBs) 
  {

    if(device_buffer_generator != nullptr or host_buffer_generator != nullptr ) {
      app_error()<<" Buffer generators already initialized in initialize_buffer_generator." 
                 <<std::endl;
      APP_ABORT(""); 
    }
//#if defined(ENABLE_ACCELERATOR)
#if defined(ENABLE_CUDA)
    device_buffer_generator = new
          BufferHandlerGenerator< device_memory_resource, device_allocator<char> >(
                      device_memory_resource{},
                      size_in_MBs*std::size_t(1024*1024),
                      device_allocator<char>{} );
    host_buffer_generator = new
          BufferHandlerGenerator< host_memory_resource, std::allocator<char> >(
                      host_memory_resource{},
                      size_in_MBs*std::size_t(1024*1024),
                      std::allocator<char>{} );
#else
    host_buffer_generator = new
          BufferHandlerGenerator< host_memory_resource, std::allocator<char> >(
                      host_memory_resource{},
                      size_in_MBs*std::size_t(1024*1024),
                      std::allocator<char>{} );
    device_buffer_generator = host_buffer_generator;
#endif
  }

  void update_buffer_generators() {
    if(device_buffer_generator == nullptr or host_buffer_generator == nullptr ) {
      app_error()<<" Uninitialized buffer generators in update_buffer_generator."
                 <<std::endl;
      APP_ABORT("");
    }    
    device_buffer_generator->update();    
    host_buffer_generator->update();    
  }  

  void destroy_buffer_generators() {
    if(host_buffer_generator != nullptr) delete host_buffer_generator;
//#if defined(ENABLE_ACCELERATOR)
#if defined(ENABLE_CUDA)
    if(device_buffer_generator != nullptr) delete device_buffer_generator;
#endif
  }
*/

}
}

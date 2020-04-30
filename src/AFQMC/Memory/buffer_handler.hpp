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

#ifndef AFQMC_BUFFER_HANDLER_HPP
#define AFQMC_BUFFER_HANDLER_HPP

#include "AFQMC/config.h"
#include <cstddef>

// new allocators
#include "multi/memory/fallback.hpp"
#include "multi/memory/allocator.hpp"
#include "multi/memory/stack.hpp"

namespace qmcplusplus
{
namespace afqmc
{

  // should I take Align from Alloc???
  template<class Alloc, class Constructor = std::allocator<char>, 
           std::size_t Align = alignof(std::max_align_t)>
  class BufferHandler
  {  

    private:

      using Ptr = typename Alloc<char>::pointer;
      using stack_mr = boost::multi::memory::stack<Ptr,Align>;
      using fallback = boost::multi::memory::fallback<stack_mr,Alloc<char>>;
      using _alloc_ = boost::multi::memory::allocator;

      Alloc base_alloc; 
      std::size_t _size=0;
      Ptr _start=nullptr;
      fallback mr_;  
      Constructor constr_;

    public:

      template<class T> 
        using allocator = _alloc_<T, fallback, 
                              typename Constructor::template rebind<T>::other>; 

      BufferHandler( Alloc& a, int initial_size=0 , Constructor const& c = {}): 
            base_alloc(a), _size(initial_size),
            _start(base_alloc.allocate(_size)), 
            mr_({_start,_size},std::addressof(base_alloc)),
            constr_(c) 
      {}

      template<class T>
      allocator<T> get_allocator(){
        return allocator<T>{mr_};
      }

  }

}
}

#endif

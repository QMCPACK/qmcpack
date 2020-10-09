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

#include <cstddef>

#include "AFQMC/config.h"
#include "mpi3/shared_communicator.hpp"

// new allocators
#include "multi/memory/fallback.hpp"
#include "multi/memory/allocator.hpp"
#include "multi/memory/stack.hpp"

namespace qmcplusplus
{
namespace afqmc
{
// should I take Align from Alloc???
template<class MemoryResource, class Constructor = std::allocator<char>, std::size_t Align = alignof(std::max_align_t)>
class BufferAllocatorGenerator
{
private:
  using memory_type = MemoryResource;
  memory_type base_mr;
  using base_element = char; // switch to byte in c++17
  using raw_pointer  = decltype(std::declval<memory_type&>().allocate(0, 0));
  using pointer      = typename std::pointer_traits<decltype(
      std::declval<memory_type&>().allocate(0, 0))>::template rebind<base_element>;
  using stack_mr     = boost::multi::memory::stack<pointer, Align>;
  using fallback     = boost::multi::memory::fallback<stack_mr, memory_type>;

  long _size     = 0;
  pointer _start = nullptr;
  fallback mr_;
  Constructor constr_;

public:
  template<class T>
  using allocator = boost::multi::memory::allocator<T, fallback, typename Constructor::template rebind<T>::other>;

  BufferAllocatorGenerator(MemoryResource const& a, long initial_size = 0, Constructor const& c = {})
      : base_mr(a),
        _size(initial_size),
        _start(static_cast<pointer>(base_mr.allocate(_size, Align))),
        mr_({_start, _size}, std::addressof(base_mr)),
        constr_(c)
  {}

  ~BufferAllocatorGenerator()
  {
    if (_start != nullptr)
      base_mr.deallocate(static_cast<raw_pointer>(_start), _size);
  }

  // sets _size to mr_.max_needed() and creates a new mr_ with capacity _size.
  // should check that stack allocator is empty
  void update() { resize(mr_.max_needed()); }

  // sets size of buffer to "at least" new_size elements
  void resize(long new_size)
  {
    if (new_size > _size)
    {
      app_log() << "\n********************************************************* \n"
                << "   Resizing memory buffer to: " << new_size / 1024.0 / 1024.0 << " MBs. \n"
                << "********************************************************* \n"
                << std::endl;
      //          mr_.reset();
      if (_size > 0)
        base_mr.deallocate(static_cast<raw_pointer>(_start), _size);
      _size  = new_size + 1024;
      _start = static_cast<pointer>(base_mr.allocate(_size, Align));
      // useful to set to zero in GPUs
      using std::fill_n;
      fill_n(_start, _size, base_element(0));
      mr_ = fallback{{_start, _size}, std::addressof(base_mr)};
    }
  }

  template<class T>
  allocator<T> get_allocator()
  {
    return allocator<T>{std::addressof(mr_), constr_};
  }
};

} // namespace afqmc
} // namespace qmcplusplus

#endif

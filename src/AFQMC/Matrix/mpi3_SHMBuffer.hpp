////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////


#ifndef  AFQMC_MPI3_SHMBUFFER_HPP 
#define  AFQMC_MPI3_SHMBUFFER_HPP 

#include "mpi3/shared_window.hpp"
#include "mpi3/shared_communicator.hpp"
#include "mpi3/shm/allocator.hpp"
#include "multi/array.hpp"
#include <memory>

namespace qmcplusplus
{

namespace afqmc
{

template<typename T>
class mpi3_SHMBuffer
{
  using communicator = boost::mpi3::shared_communicator;
  using shared_allocator = boost::mpi3::shm::allocator<T>;
  using shared_array = boost::multi::array<T,1,shared_allocator>;

  public:

    mpi3_SHMBuffer(communicator& comm_, size_t n=0):
        comm(std::addressof(comm_)),array_({n},shared_allocator{comm_})
    {
    }

    mpi3_SHMBuffer<T>(const mpi3_SHMBuffer<T>& other) = delete;
    mpi3_SHMBuffer<T>& operator=(const mpi3_SHMBuffer<T>& other) = delete;

    mpi3_SHMBuffer<T>(mpi3_SHMBuffer<T>&& other) = default; 
    mpi3_SHMBuffer<T>& operator=(mpi3_SHMBuffer<T>&& other) = default; 
/*
    { 
      if(this != &other) {
        comm = other.comm;
        other.comm = NULL;
        win = std::move(other.win);
        comm->barrier();
      }
      return *this;
    }
*/

    void resize(size_t n) {
      if(array_.size() != n)
        array_.reextent({n});
    }

    T* data() {return std::addressof(*array_.data());} 
    T const* data() const{return std::addressof(*array_.data());} 
    size_t size() const{return static_cast<size_t>(array_.size());} 

    communicator& getCommunicator() const{return *comm;}

  private:

    // mpi3_SHMBuffer does not own the communicator comm points to.
    // This class assumes that objects of this class live in the same scope
    // as the communicator comm points to and  that it will remain valid and equivalent. 
    communicator* comm;
    
    shared_array array_;

};

}

}

#endif

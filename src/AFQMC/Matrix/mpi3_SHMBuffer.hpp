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
#include <memory>

namespace qmcplusplus
{

namespace afqmc
{

template<typename T>
class mpi3_SHMBuffer
{
  using communicator = boost::mpi3::shared_communicator;
  using shared_window = boost::mpi3::shared_window<T>;  

  public:

    mpi3_SHMBuffer(communicator& comm_, size_t n=0):
        comm(std::addressof(comm_)),win(comm_,(comm_.root()?n:0),sizeof(T))
    {
    }

    mpi3_SHMBuffer<T>(const mpi3_SHMBuffer<T>& other) = delete;
    mpi3_SHMBuffer<T>& operator=(const mpi3_SHMBuffer<T>& other) = delete;

    mpi3_SHMBuffer<T>(mpi3_SHMBuffer<T>&& other) {
      *this = std::move(other);
    }

    mpi3_SHMBuffer<T>& operator=(mpi3_SHMBuffer<T>&& other) {
      if(this != &other) {
        comm = other.comm;
        other.comm = NULL;
        win = std::move(other.win);
        comm->barrier();
      }
      return *this;
    }

    void resize(size_t n) {
      if(size() == n) return;
      shared_window w0(*comm,(comm->root()?n:0),sizeof(T)); 
      comm->barrier();
      if(comm->rank()==0) {
        if(size() < n) 
          std::copy(this->data(),this->data()+size(),w0.base(0));
        else
          std::copy(this->data(),this->data()+static_cast<size_t>(w0.size(0)),w0.base(0));
      }   
      comm->barrier();
      win = std::move(w0);
      comm->barrier();
    }

    T* data() {return win.base(0);} 
    T const* data() const{return win.base(0);} 
    size_t size() const{return static_cast<size_t>(win.size(0));} 

    communicator& getCommunicator() const{return *comm;}

  private:

    // mpi3_SHMBuffer does not own the communicator comm points to.
    // This class assumes that objects of this class live in the same scope
    // as the communicator comm points to and  that it will remain valid and equivalent. 
    communicator* comm;
    
    shared_window win;

};

}

}

#endif

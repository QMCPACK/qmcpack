////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////


#ifndef  AFQMC_MPI3_SHARED_MA_PROXY_HPP
#define  AFQMC_MPI3_SHARED_MA_PROXY_HPP 

#include "mpi3/shared_communicator.hpp"

namespace qmcplusplus
{

namespace afqmc
{

template<typename T>
class mpi3_shared_ma_proxy
{
  using communicator = boost::mpi3::shared_communicator;
  using shmCMatrix = boost::multi::array<T,2,shared_allocator<T>>;
  using size_type = std::size_t;

  public:

    using value_type = T;
    using ma_type = boost::multi::array_ref<T,2>;

    mpi3_shared_ma_proxy(communicator& comm_,
                         std::array<size_type, 2> s_, 
                         std::array<size_type, 2> gs_={0,0},
                         std::array<size_type, 2> o_={0,0}):   
       global_shape_(gs_),
       offset_(o_),
       M({s_[0],s_[1]},shared_allocator<T>{comm_}) 
    {
      if(gs_[0] == 0 || gs_[1] == 0) 
        global_shape_ = s_;
    }

    mpi3_shared_ma_proxy(mpi3_shared_ma_proxy const& other) = delete;
    mpi3_shared_ma_proxy& operator=(mpi3_shared_ma_proxy const& other) = delete;

    mpi3_shared_ma_proxy(mpi3_shared_ma_proxy && other) = default;
    mpi3_shared_ma_proxy& operator=(mpi3_shared_ma_proxy && other) = default;

    T* origin() {return std::addressof(*M.origin());} 
    T const* origin() const{return std::addressof(*M.origin());} 
    T* data() {return std::addressof(*M.data());} 
    T const* data() const{return std::addressof(*M.data());} 
    size_type num_elements() const{return M.num_elements(); } 
    auto shape() const{return M.shape();} 
    auto strides() const{return M.strides();} 
    std::array<size_type, 2> global_offset() const{return offset_;} 
    std::array<size_type, 2> global_shape() const{return global_shape_;} 

    boost::multi::array_ref<T,2> get() { 
        return boost::multi::array_ref<T,2>(this->origin(),M.extensions()); 
    }

  private:

    std::array<size_type, 2> global_shape_;
    std::array<size_type, 2> offset_;
    shmCMatrix M;

};

}

}

#endif

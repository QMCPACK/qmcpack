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
#include "AFQMC/Matrix/mpi3_SHMBuffer.hpp"

namespace qmcplusplus
{

namespace afqmc
{

template<typename T>
class mpi3_shared_ma_proxy
{
  using communicator = boost::mpi3::shared_communicator;
  using buffer = mpi3_SHMBuffer<T>; 
  using size_type = std::size_t;

  public:

    using value_type = T;
    using ma_type = boost::multi_array_ref<T,2>;

    mpi3_shared_ma_proxy(communicator& comm_,
                         std::array<size_type, 2> s_, 
                         std::array<size_type, 2> gs_={0,0},
                         std::array<size_type, 2> o_={0,0}):   
       shape_(s_),
       global_shape_(gs_),
       offset_(o_),
       ptr(std::make_unique<buffer>(comm_,s_[0]*s_[1])) 
    {
      if(gs_[0] == 0 || gs_[1] == 0) 
        global_shape_ = shape_;
    }

    mpi3_shared_ma_proxy(mpi3_shared_ma_proxy const& other) = delete;
    mpi3_shared_ma_proxy& operator=(mpi3_shared_ma_proxy const& other) = delete;

    mpi3_shared_ma_proxy(mpi3_shared_ma_proxy && other) = default;
    mpi3_shared_ma_proxy& operator=(mpi3_shared_ma_proxy && other) = default;

    void reshape(std::array<size_type, 2> s_, 
                 std::array<size_type, 2> gs_={0,0},
                 std::array<size_type, 2> o_={0,0}) {
      shape_=s_;
      if(gs_[0] == 0 || gs_[1] == 0) 
        global_shape_ = shape_;
      else
        global_shape_=gs_;
      offset_=o_;  
      ptr->resize(s_[0]*s_[1]);
    }

    T* origin() {return ptr->data();} 
    T const* origin() const{return ptr->data();} 
    T* data() {return ptr->data();} 
    T const* data() const{return ptr->data();} 
    size_type num_elements() const{return shape_[0]*shape_[1];} 
    std::array<size_type, 2> shape() const{return shape_;} 
    std::array<size_type, 2> strides() const{return std::array<size_type, 2>{shape_[1],1};} 
    std::array<size_type, 2> offset() const{return offset_;} 
    std::array<size_type, 2> global_shape() const{return global_shape_;} 

    boost::multi_array_ref<T,2> get() { 
        return boost::multi_array_ref<T,2>(ptr->data(),extents[shape_[0]][shape_[1]]); 
    }

  private:

    std::array<size_type, 2> shape_;
    std::array<size_type, 2> global_shape_;
    std::array<size_type, 2> offset_;
    std::unique_ptr<buffer> ptr;

};

}

}

#endif

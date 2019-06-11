//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_WALKERAVERAGED_OBSERVABLES_HPP
#define QMCPLUSPLUS_AFQMC_WALKERAVERAGED_OBSERVABLES_HPP

#include <iostream>
#include <vector>
#include <string>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class WalkerAveragedObservables: public AFQMCInfo
{

  // allocators
  using Allocator = device_allocator<ComplexType>;
  using Allocator_shared = localTG_allocator<ComplexType>;

  // type defs
  using pointer = typename Allocator::pointer;
  using const_pointer = typename Allocator::const_pointer;
  using pointer_shared = typename Allocator_shared::pointer;
  using const_pointer_shared = typename Allocator_shared::const_pointer;

  using CVector = boost::multi::array<ComplexType,1,Allocator>;
  using CMatrix = boost::multi::array<ComplexType,2,Allocator>;
  using CTensor = boost::multi::array<ComplexType,3,Allocator>;
  using CVector_ref = boost::multi::array_ref<ComplexType,1,pointer>;
  using CMatrix_ref = boost::multi::array_ref<ComplexType,2,pointer>;
  using CMatrix_cref = boost::multi::array_ref<const ComplexType,2,const_pointer>;
  using CTensor_ref = boost::multi::array_ref<ComplexType,3,pointer>;
  using CTensor_cref = boost::multi::array_ref<const ComplexType,3,const_pointer>;
  using shmCVector = boost::multi::array<ComplexType,1,Allocator_shared>;
  using shared_mutex = boost::mpi3::shm::mutex;

  using stdCVector = boost::multi::array<ComplexType,1>;
  using stdCMatrix = boost::multi::array<ComplexType,2>;
  using stdCTensor = boost::multi::array<ComplexType,3>;
  using mpi3CVector = boost::multi::array<ComplexType,1,shared_allocator<ComplexType>>;
  using mpi3CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using mpi3CTensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;

  using stdCMatrix_ref = boost::multi::array_ref<ComplexType,2>;


  public: 




  private:

    Wavefunction& ref_wfn;   
    
    SlaterDetOperations* SDet;

     


}; 

}

}

#endif

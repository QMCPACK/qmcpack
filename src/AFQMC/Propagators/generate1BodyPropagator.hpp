//////////////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_GENERATEPROPAGATOR_HPP
#define QMCPLUSPLUS_AFQMC_GENERATEPROPAGATOR_HPP

#include<cstdlib>
#include<algorithm>
#include<complex>
#include<iostream>
#include<vector>
#include<numeric>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "Configuration.h"
#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/Matrix/csr_matrix_construct.hpp"

namespace qmcplusplus
{
namespace afqmc
{

  // Input H1(i,j) = h(i,j) + sum_n vMF(n)*CholMat(i,j,n) + vn0(i,j)  
  // Output: sparse 1 body propagator = exp(-0.5*dt*H1) 
  template<class MultiArray2D>  
  P1Type generate1BodyPropagator(TaskGroup_& TG, RealType cut, RealType dt, MultiArray2D const& H1)
  {
    assert(H1.dimensionality==2);
    assert(H1.shape()[0] == H1.shape()[1]);
    assert(H1.strides()[1] == 1);
    int NMO = H1.shape()[0];
    if(TG.TG_local().root()) {   
      boost::multi_array<ComplexType,2> v(extents[NMO][NMO]);
      std::fill_n(v.origin(),v.num_elements(),ComplexType(0));  

      for(int i=0; i<NMO; i++) 
        ma::axpy(-0.5*dt,H1[i],v[i]);

      boost::multi_array<ComplexType,2> P = ma::exp(v);

      return csr::shm::construct_csr_matrix_single_input<P1Type>(P,cut,'N',TG.TG_local());
    } else {
      boost::multi_array<ComplexType,2> P(extents[1][1]);
      return csr::shm::construct_csr_matrix_single_input<P1Type>(P,cut,'N',TG.TG_local());
    }    
  }

}
}
#endif

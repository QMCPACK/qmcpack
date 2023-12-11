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

#include <cstdlib>
#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>
#include <numeric>
#if defined(USE_MPI)
#include <mpi.h>
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
template<class P_Type, class MultiArray2D>
P_Type generate1BodyPropagator(TaskGroup_& TG,
                               RealType cut,
                               RealType dt,
                               MultiArray2D const& H1,
                               bool printP1eV = false)
{
  assert(H1.dimensionality == 2);
  assert(std::get<0>(H1.sizes()) == std::get<1>(H1.sizes()));
  assert(H1.stride(1) == 1);
  int NMO = H1.size();
  if (TG.TG_local().root())
  {
    boost::multi::array<ComplexType, 2> v({NMO, NMO});
    fill_n(v.origin(), v.num_elements(), ComplexType(0));

    // running on host regardless
    boost::multi::array<ComplexType, 2> h1_(H1.extensions());
    std::copy_n(to_address(H1.origin()), H1.num_elements(), h1_.origin());

    for (int i = 0; i < NMO; i++)
      ma::axpy(-0.5 * dt, h1_[i], v[i]);

    boost::multi::array<ComplexType, 2> P = ma::exp(v, printP1eV);

    // need a version of this that works with gpu_ptr!!!
    return csr::shm::construct_csr_matrix_single_input<P_Type>(P, cut, 'N', TG.TG_local());
  }
  else
  {
    boost::multi::array<ComplexType, 2> P({1, 1});
    return csr::shm::construct_csr_matrix_single_input<P_Type>(P, cut, 'N', TG.TG_local());
  }
}

// Input H1(i,j) = h(i,j) + sum_n vMF(n)*CholMat(i,j,n) + vn0(i,j)
// Output: sparse 1 body propagator = exp(-0.5*dt*H1)
template<class P_Type, class MultiArray2D, class MultiArray2DB>
P_Type generate1BodyPropagator(TaskGroup_& TG,
                               RealType cut,
                               RealType dt,
                               MultiArray2D const& H1,
                               MultiArray2DB const& H1ext,
                               bool printP1eV = false)
{
  assert(H1.dimensionality == 2);
  assert(std::get<0>(H1.sizes()) == std::get<1>(H1.sizes()));
  assert(H1.stride(1) == 1);
  assert(H1ext.dimensionality == 2);
  assert(std::get<0>(H1ext.sizes()) == std::get<1>(H1ext.sizes()));
  assert(H1ext.stride(1) == 1);
  assert(std::get<0>(H1.sizes()) == std::get<1>(H1ext.sizes()));
  int NMO = H1.size();
  if (TG.TG_local().root())
  {
    //      boost::multi::array<ComplexType,2> v({NMO,NMO});
    //      fill_n(v.origin(),v.num_elements(),ComplexType(0));

    // running on host regardless
    boost::multi::array<ComplexType, 2> h1_(H1);
    //boost::multi::array<ComplexType,2> h1ext_(H1ext);
    boost::multi::array<ComplexType, 2> h1ext_({NMO, NMO});
    //copy_n(H1ext.origin(),NMO*NMO,h1ext_.origin());
    h1ext_ = H1ext;


    ma::add(ComplexType(-0.5 * dt), h1_, ComplexType(-0.5 * dt), h1ext_, h1_);

    //      for(int i=0; i<NMO; i++)
    //        ma::axpy(-0.5*dt,h1_[i],v[i]);

    boost::multi::array<ComplexType, 2> P = ma::exp(h1_, printP1eV);

    // need a version of this that works with gpu_ptr!!!
    return csr::shm::construct_csr_matrix_single_input<P_Type>(P, cut, 'N', TG.TG_local());
  }
  else
  {
    boost::multi::array<ComplexType, 2> P({1, 1});
    return csr::shm::construct_csr_matrix_single_input<P_Type>(P, cut, 'N', TG.TG_local());
  }
}

} // namespace afqmc
} // namespace qmcplusplus
#endif

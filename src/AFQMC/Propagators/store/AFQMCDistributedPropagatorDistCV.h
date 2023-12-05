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

#ifndef QMCPLUSPLUS_AFQMC_DISTRIBUTEDPROPAGATORDISTCV_H
#define QMCPLUSPLUS_AFQMC_DISTRIBUTEDPROPAGATORDISTCV_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "hdf/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "mpi3/shm/mutex.hpp"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"

#include "AFQMC/Propagators/AFQMCBasePropagator.h"

namespace qmcplusplus
{
namespace afqmc
{
/*
 * AFQMC propagator using distributed Cholesky matrices over nodes in the TG. 
 * Shared memory parallelization is used for on-node concurrency.
 *   - Specialized algorithm for case when vbias doesn't need to be reduced over 
 *     nodes in a TG.
 */
class AFQMCDistributedPropagatorDistCV : public AFQMCBasePropagator
{
  using base            = AFQMCBasePropagator;
  using stdCVector      = ComplexVector<std::allocator<ComplexType>>;
  using stdCVector_ref  = ComplexVector_ref<ComplexType*>;
  using stdCMatrix_ref  = ComplexMatrix_ref<ComplexType*>;
  using stdC3Tensor_ref = boost::multi::array_ref<ComplexType, 3>;

public:
  AFQMCDistributedPropagatorDistCV(AFQMCInfo& info,
                                   xmlNodePtr cur,
                                   afqmc::TaskGroup_& tg_,
                                   Wavefunction& wfn_,
                                   CMatrix&& h1_,
                                   CVector&& vmf_,
                                   RandomBase<RealType>& r)
      : base(info, cur, tg_, wfn_, std::move(h1_), std::move(vmf_), r),
        req_Gsend(MPI_REQUEST_NULL),
        req_Grecv(MPI_REQUEST_NULL),
        req_vsend(MPI_REQUEST_NULL),
        req_vrecv(MPI_REQUEST_NULL)
  {
    assert(TG.getNGroupsPerTG() > 1);
  }

  ~AFQMCDistributedPropagatorDistCV()
  {
    if (req_Grecv != MPI_REQUEST_NULL)
      MPI_Request_free(&req_Grecv);
    if (req_Gsend != MPI_REQUEST_NULL)
      MPI_Request_free(&req_Gsend);
    if (req_vrecv != MPI_REQUEST_NULL)
      MPI_Request_free(&req_vrecv);
    if (req_vsend != MPI_REQUEST_NULL)
      MPI_Request_free(&req_vsend);
  }

  AFQMCDistributedPropagatorDistCV(AFQMCDistributedPropagatorDistCV const& other) = delete;
  AFQMCDistributedPropagatorDistCV& operator=(AFQMCDistributedPropagatorDistCV const& other) = delete;
  AFQMCDistributedPropagatorDistCV(AFQMCDistributedPropagatorDistCV&& other)                 = default;
  AFQMCDistributedPropagatorDistCV& operator=(AFQMCDistributedPropagatorDistCV&& other) = default;

  template<class WlkSet>
  void Propagate(int steps, WlkSet& wset, RealType E1, RealType dt, int fix_bias = 1)
  {
    int nblk   = steps / fix_bias;
    int nextra = steps % fix_bias;
    for (int i = 0; i < nblk; i++)
      step(fix_bias, wset, E1, dt);
    if (nextra > 0)
      step(nextra, wset, E1, dt);
    TG.local_barrier();
  }

protected:
  MPI_Request req_Gsend, req_Grecv;
  MPI_Request req_vsend, req_vrecv;

#ifdef ENABLE_CUDA
  stdCVector comm_buffer;
#else
  sharedCVector comm_buffer;
#endif

  template<class WlkSet>
  void step(int steps, WlkSet& wset, RealType E1, RealType dt);
};

} // namespace afqmc

} // namespace qmcplusplus

#include "AFQMC/Propagators/AFQMCDistributedPropagatorDistCV.icc"

#endif

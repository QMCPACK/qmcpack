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
#include "AFQMC/Memory/buffer_managers.h"

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
  using base = AFQMCBasePropagator;

public:
  AFQMCDistributedPropagatorDistCV(AFQMCInfo& info,
                                   xmlNodePtr cur,
                                   afqmc::TaskGroup_& tg_,
                                   Wavefunction& wfn_,
                                   stdCMatrix&& h1_,
                                   CVector&& vmf_,
                                   RandomBase<RealType>& r)
      : base(info, cur, tg_, wfn_, std::move(h1_), std::move(vmf_), r),
        bpX(iextensions<1u>{1}, shared_allocator<ComplexType>{TG.TG_local()}),
        req_Gsend(MPI_REQUEST_NULL),
        req_Grecv(MPI_REQUEST_NULL),
        req_vsend(MPI_REQUEST_NULL),
        req_vrecv(MPI_REQUEST_NULL),
        req_Xsend(MPI_REQUEST_NULL),
        req_Xrecv(MPI_REQUEST_NULL),
        req_X2send(MPI_REQUEST_NULL),
        req_X2recv(MPI_REQUEST_NULL),
        req_bpvsend(MPI_REQUEST_NULL),
        req_bpvrecv(MPI_REQUEST_NULL)
  {
    assert(TG.getNGroupsPerTG() > 1);

    std::string str("no");
    ParameterSet m_param;
    m_param.add(str, "low_memory");
    m_param.put(cur);

    std::transform(str.begin(), str.end(), str.begin(), (int (*)(int))tolower);
    if (str == "yes" || str == "true")
      low_memory_step = true;

    if (low_memory_step)
      app_log() << " Using low memory distributed propagation. \n";
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
    if (req_X2recv != MPI_REQUEST_NULL)
      MPI_Request_free(&req_X2recv);
    if (req_X2send != MPI_REQUEST_NULL)
      MPI_Request_free(&req_X2send);
    if (req_Xrecv != MPI_REQUEST_NULL)
      MPI_Request_free(&req_Xrecv);
    if (req_Xsend != MPI_REQUEST_NULL)
      MPI_Request_free(&req_Xsend);
    if (req_bpvrecv != MPI_REQUEST_NULL)
      MPI_Request_free(&req_bpvrecv);
    if (req_bpvsend != MPI_REQUEST_NULL)
      MPI_Request_free(&req_bpvsend);
  }

  AFQMCDistributedPropagatorDistCV(AFQMCDistributedPropagatorDistCV const& other) = delete;
  AFQMCDistributedPropagatorDistCV& operator=(AFQMCDistributedPropagatorDistCV const& other) = delete;
  AFQMCDistributedPropagatorDistCV(AFQMCDistributedPropagatorDistCV&& other)                 = default;
  AFQMCDistributedPropagatorDistCV& operator=(AFQMCDistributedPropagatorDistCV&& other) = delete;

  template<class WlkSet>
  void Propagate(int steps, WlkSet& wset, RealType E1, RealType dt, int fix_bias = 1)
  {
    int nblk   = steps / fix_bias;
    int nextra = steps % fix_bias;
    if (low_memory_step)
    {
      for (int i = 0; i < nblk; i++)
      {
        step_collective(fix_bias, wset, E1, dt);
        update_memory_managers();
      }
      if (nextra > 0)
        step_collective(nextra, wset, E1, dt);
    }
    else
    {
      for (int i = 0; i < nblk; i++)
      {
        step(fix_bias, wset, E1, dt);
        update_memory_managers();
      }
      if (nextra > 0)
        step(nextra, wset, E1, dt);
    }
    TG.local_barrier();
  }

  template<class WlkSet, class CTens, class CMat>
  void BackPropagate(int steps, int nStabalize, WlkSet& wset, CTens&& Refs, CMat&& detR);

protected:
  mpi3SPCVector bpX;
  std::vector<int> bpx_counts, bpx_displ;

  bool buffer_reallocated    = false;
  bool buffer_reallocated_bp = false;
  bool low_memory_step       = false;

  MPI_Request req_Gsend, req_Grecv;
  MPI_Request req_vsend, req_vrecv;

  MPI_Request req_Xsend, req_Xrecv;
  MPI_Request req_X2send, req_X2recv;
  MPI_Request req_bpvsend, req_bpvrecv;

  template<class WlkSet>
  void step(int steps, WlkSet& wset, RealType E1, RealType dt);

  template<class WlkSet>
  void step_collective(int steps, WlkSet& wset, RealType E1, RealType dt);
};

} // namespace afqmc

} // namespace qmcplusplus

#include "AFQMC/Propagators/AFQMCDistributedPropagatorDistCV.icc"

#endif

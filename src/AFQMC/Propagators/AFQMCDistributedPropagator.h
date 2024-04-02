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

#ifndef QMCPLUSPLUS_AFQMC_DISTRIBUTEDPROPAGATOR_H
#define QMCPLUSPLUS_AFQMC_DISTRIBUTEDPROPAGATOR_H

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
 *  - General case. Both vbias and vHS are assumed to need reduction over all nodes
 *    in TG.
 */
class AFQMCDistributedPropagator : public AFQMCBasePropagator
{
  using base = AFQMCBasePropagator;

public:
  AFQMCDistributedPropagator(AFQMCInfo& info,
                             xmlNodePtr cur,
                             afqmc::TaskGroup_& tg_,
                             Wavefunction& wfn_,
                             stdCMatrix&& h1_,
                             CVector&& vmf_,
                             RandomBase<RealType>& r)
      : base(info, cur, tg_, wfn_, std::move(h1_), std::move(vmf_), r),
        core_comm(tg_.TG().split(tg_.getLocalTGRank(), tg_.TG().rank()))
  //            ,core_comm()
  {
    //      core_comm = std::move(tg_.TG().split(tg_.getLocalTGRank()));
    assert(TG.getNGroupsPerTG() > 1);
  }

  ~AFQMCDistributedPropagator() {}

  AFQMCDistributedPropagator(AFQMCDistributedPropagator const& other) = delete;
  AFQMCDistributedPropagator& operator=(AFQMCDistributedPropagator const& other) = delete;
  //AFQMCDistributedPropagator(AFQMCDistributedPropagator&& other) = default;
  AFQMCDistributedPropagator(AFQMCDistributedPropagator&& other) : base(std::move(other)), core_comm()
  {
    // move constructor for communicator seems broken
    core_comm = TG.TG().split(TG.getLocalTGRank(), TG.TG().rank());
  }
  AFQMCDistributedPropagator& operator=(AFQMCDistributedPropagator&& other) = delete;

  template<class WlkSet>
  void Propagate(int steps, WlkSet& wset, RealType E1, RealType dt, int fix_bias = 1)
  {
    int nblk   = steps / fix_bias;
    int nextra = steps % fix_bias;
    for (int i = 0; i < nblk; i++)
    {
      step(fix_bias, wset, E1, dt);
      update_memory_managers();
    }
    if (nextra > 0)
      step(nextra, wset, E1, dt);
    TG.local_barrier();
  }

  template<class WlkSet, class CTens, class CMat>
  void BackPropagate(int steps, int nStabalize, WlkSet& wset, CTens&& Refs, CMat&& detR)
  {
    APP_ABORT(" Error: Finish BackPropagate.\n");
  }

protected:
  // new communicator over similar cores in a TG
  // every core communicates a segment to increase effective bandwidth
  boost::mpi3::communicator core_comm;

  // additional dimension for temporary computation
  C3Tensor MFfactor;
  C3Tensor hybrid_weight;

  template<class WlkSet>
  void step(int steps, WlkSet& wset, RealType E1, RealType dt);
};

} // namespace afqmc

} // namespace qmcplusplus

#include "AFQMC/Propagators/AFQMCDistributedPropagator.icc"

#endif

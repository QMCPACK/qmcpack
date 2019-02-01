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

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "mpi3/shm/mutex.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"

#include "AFQMC/Propagators/AFQMCSharedPropagator.h"

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
class AFQMCDistributedPropagator: public AFQMCSharedPropagator
{
  using base = AFQMCSharedPropagator;
  public:

    AFQMCDistributedPropagator(AFQMCInfo& info, xmlNodePtr cur, afqmc::TaskGroup_& tg_, 
                          Wavefunction& wfn_, CMatrix&& h1_, CVector&& vmf_, 
                          RandomGenerator_t* r): 
            AFQMCSharedPropagator(info,cur,tg_,wfn_,std::move(h1_),std::move(vmf_),r)
            ,core_comm(tg_.TG().split(tg_.getLocalTGRank(),tg_.TG().rank()))
//            ,core_comm()
    {
//      core_comm = std::move(tg_.TG().split(tg_.getLocalTGRank()));
      assert(TG.getNNodesPerTG() > 1);
    }

    ~AFQMCDistributedPropagator() {}

    AFQMCDistributedPropagator(AFQMCDistributedPropagator const& other) = delete;
    AFQMCDistributedPropagator& operator=(AFQMCDistributedPropagator const& other) = delete;
    //AFQMCDistributedPropagator(AFQMCDistributedPropagator&& other) = default;
    AFQMCDistributedPropagator(AFQMCDistributedPropagator&& other):
                AFQMCSharedPropagator(std::move(other)),
                core_comm()
    {
      // move constructor for communicator seems broken
      core_comm = std::move(TG.TG().split(TG.getLocalTGRank(),TG.TG().rank()));  
    }
    AFQMCDistributedPropagator& operator=(AFQMCDistributedPropagator&& other) = default;

    template<class WlkSet>
    void Propagate(int steps, WlkSet& wset, RealType E1,
                   RealType dt, int fix_bias=1) {
      int nblk = steps/fix_bias;
      int nextra = steps%fix_bias;
      for(int i=0; i<nblk; i++)
        step(fix_bias,wset,E1,dt);
      if(nextra>0)
        step(nextra,wset,E1,dt);
      TG.local_barrier();
    }

  protected: 

    // new communicator over similar cores in a TG
    // every core communicates a segment to increase effective bandwidth
    boost::mpi3::communicator core_comm;

    template<class WlkSet>
    void step(int steps, WlkSet& wset, RealType E1, RealType dt);

    // additional dimension for temporary computation
    boost::multi::array<ComplexType,3> MFfactor;
    boost::multi::array<ComplexType,3> hybrid_weight;

};

}

}

#include "AFQMC/Propagators/AFQMCDistributedPropagator.icc"

#endif


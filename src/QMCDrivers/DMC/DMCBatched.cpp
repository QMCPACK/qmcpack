//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: DMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/DMC/DMCBatched.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Concurrency/Info.hpp"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"


namespace qmcplusplus
{
/** Constructor maintains proper ownership of input parameters
   */
DMCBatched::DMCBatched(QMCDriverInput&& qmcdriver_input,
                       DMCDriverInput&& input,
                       MCPopulation&& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& wf_pool,
                       Communicate* comm)
    : QMCDriverNew(std::move(qmcdriver_input), std::move(pop), psi, h, wf_pool, "DMCBatched::", comm), 
      dmcdriver_input_(input)
{
  QMCType = "DMCBatched";
}

QMCTraits::IndexType DMCBatched::calc_default_local_walkers(IndexType walkers_per_rank)
{
  checkNumCrowdsLTNumThreads();
  int num_threads(Concurrency::maxThreads<>());
  IndexType rw = walkers_per_rank; //qmcdriver_input_.get_walkers_per_rank();
  if (num_crowds_ == 0)
    num_crowds_ = std::min(num_threads, rw);
  walkers_per_crowd_      = (rw % num_crowds_) ? rw / num_crowds_ + 1 : rw / num_crowds_;

  IndexType local_walkers = walkers_per_crowd_ * num_crowds_;
  population_.set_num_local_walkers(local_walkers);
  population_.set_num_global_walkers(local_walkers * population_.get_num_ranks());
  if (rw != qmcdriver_input_.get_walkers_per_rank())
    app_warning() << "DMCBatched driver has adjusted walkers per rank to: " << local_walkers << '\n';

  app_log() << "DMCBatched walkers per crowd " << walkers_per_crowd_ << std::endl;
  return local_walkers;
}

bool DMCBatched::run()
{
  LoopTimer dmc_loop;

  return false;
}

} // namespace qmcplusplus

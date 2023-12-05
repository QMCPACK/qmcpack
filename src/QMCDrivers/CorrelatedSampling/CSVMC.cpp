//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "CSVMC.h"
#include "QMCDrivers/CorrelatedSampling/CSVMCUpdateAll.h"
#include "QMCDrivers/CorrelatedSampling/CSVMCUpdatePbyP.h"
#include "Estimators/CSEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"
#include "RandomNumberControl.h"
#include "Concurrency/OpenMP.h"
#include "Message/CommOperators.h"
#include "Utilities/FairDivide.h"
#include "Utilities/qmc_common.h"
//#define ENABLE_VMC_OMP_MASTER
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif

namespace qmcplusplus
{
/// Constructor.
CSVMC::CSVMC(const ProjectData& project_data,
             MCWalkerConfiguration& w,
             TrialWaveFunction& psi,
             QMCHamiltonian& h,
             Communicate* comm)
    : QMCDriver(project_data, w, psi, h, comm, "CSVMC"), UseDrift("yes"), multiEstimator(0), Mover(0)
{
  RootName = "csvmc";
  m_param.add(UseDrift, "useDrift");
  m_param.add(UseDrift, "usedrift");
  m_param.add(UseDrift, "use_drift");
  equilBlocks = -1;
  m_param.add(equilBlocks, "equilBlocks");
  qmc_driver_mode.set(QMC_MULTIPLE, 1);
}

/** allocate internal data here before run() is called
 * @author SIMONE
 *
 * See QMCDriver::process
 */
bool CSVMC::put(xmlNodePtr q)
{
  int target_min = -1;
  ParameterSet p;
  p.add(target_min, "minimumtargetwalkers");
  p.add(target_min, "minimumsamples");
  p.put(q);

  app_log() << "\n<vmc function=\"put\">"
            << "\n  qmc_counter=" << qmc_common.qmc_counter << "  my_counter=" << MyCounter << std::endl;
  if (qmc_common.qmc_counter && MyCounter)
  {
    nSteps               = prevSteps;
    nStepsBetweenSamples = prevStepsBetweenSamples;
  }
  else
  {
    int nw = W.getActiveWalkers();
    //compute samples and overwrite steps for the given samples
    int Nthreads = omp_get_max_threads();
    int Nprocs   = myComm->size();
    //target samples set by samples or samplesperthread/dmcwalkersperthread
    nTargetPopulation = std::max(nTargetPopulation, nSamplesPerThread * Nprocs * Nthreads);
    nTargetSamples    = static_cast<int>(std::ceil(nTargetPopulation));

    if (nTargetSamples)
    {
      int nwtot      = nw * Nprocs; //total number of walkers used by this qmcsection
      nTargetSamples = std::max(nwtot, nTargetSamples);
      if (target_min > 0)
      {
        nTargetSamples    = std::max(nTargetSamples, target_min);
        nTargetPopulation = std::max(nTargetPopulation, static_cast<RealType>(target_min));
      }
      nTargetSamples = ((nTargetSamples + nwtot - 1) / nwtot) *
          nwtot; // nTargetSamples are always multiples of total number of walkers
      nSamplesPerThread = nTargetSamples / Nprocs / Nthreads;
      int ns_target     = nTargetSamples * nStepsBetweenSamples; //total samples to generate
      int ns_per_step   = Nprocs * nw;                           //total samples per step
      nSteps            = std::max(nSteps, (ns_target / ns_per_step + nBlocks - 1) / nBlocks);
      Period4WalkerDump = nStepsBetweenSamples = (ns_per_step * nSteps * nBlocks) / nTargetSamples;
    }
    else
    {
      Period4WalkerDump = nStepsBetweenSamples = (nBlocks + 1) * nSteps; //some positive number, not used
      nSamplesPerThread                        = 0;
    }
  }
  prevSteps               = nSteps;
  prevStepsBetweenSamples = nStepsBetweenSamples;

  app_log() << "  time step      = " << Tau << std::endl;
  app_log() << "  blocks         = " << nBlocks << std::endl;
  app_log() << "  steps          = " << nSteps << std::endl;
  app_log() << "  substeps       = " << nSubSteps << std::endl;
  app_log() << "  current        = " << CurrentStep << std::endl;
  app_log() << "  target samples = " << nTargetPopulation << std::endl;
  app_log() << "  walkers/mpi    = " << W.getActiveWalkers() << std::endl << std::endl;
  app_log() << "  stepsbetweensamples = " << nStepsBetweenSamples << std::endl;

  m_param.get(app_log());

  if (DumpConfig)
  {
    app_log() << "  DumpConfig==true Configurations are dumped to config.h5 with a period of " << Period4CheckPoint
              << " blocks" << std::endl;
  }
  else
  {
    app_log() << "  DumpConfig==false Nothing (configurations, state) will be saved." << std::endl;
  }

  if (Period4WalkerDump > 0)
    app_log() << "  Walker Samples are dumped every " << Period4WalkerDump << " steps." << std::endl;

  app_log() << "</vmc>" << std::endl;
  app_log().flush();

  return true;
}

/** Run the CSVMC algorithm.
 *
 * Similar to VMC::run
 */
bool CSVMC::run()
{
  resetRun();
  //start the main estimator
  Estimators->start(nBlocks);
  for (int ip = 0; ip < NumThreads; ++ip)
    CSMovers[ip]->startRun(nBlocks, false);

#if !defined(REMOVE_TRACEMANAGER)
  Traces->startRun(nBlocks, traceClones);
#endif
  const bool has_collectables = W.Collectables.size();
  for (int block = 0; block < nBlocks; ++block)
  {
#pragma omp parallel
    {
      int ip                 = omp_get_thread_num();
      IndexType updatePeriod = (qmc_driver_mode[QMC_UPDATE_MODE]) ? Period4CheckProperties : 0;
      //assign the iterators and resuse them
      MCWalkerConfiguration::iterator wit(W.begin() + wPerRank[ip]), wit_end(W.begin() + wPerRank[ip + 1]);
      CSMovers[ip]->startBlock(nSteps);
      int now_loc    = CurrentStep;
      RealType cnorm = 1.0 / static_cast<RealType>(wPerRank[ip + 1] - wPerRank[ip]);
      for (int step = 0; step < nSteps; ++step)
      {
        CSMovers[ip]->set_step(now_loc);
        //collectables are reset, it is accumulated while advancing walkers
        wClones[ip]->resetCollectables();
        CSMovers[ip]->advanceWalkers(wit, wit_end, false);
        if (has_collectables)
          wClones[ip]->Collectables *= cnorm;
        CSMovers[ip]->accumulate(wit, wit_end);
        ++now_loc;
        if (Period4WalkerDump && now_loc % Period4WalkerDump == 0)
          wClones[ip]->saveEnsemble(wit, wit_end);
      }
      CSMovers[ip]->stopBlock(false);

    } //end-of-parallel for
    CurrentStep += nSteps;
    Estimators->stopBlock(estimatorClones);
#if !defined(REMOVE_TRACEMANAGER)
    Traces->write_buffers(traceClones, block);
#endif
    recordBlock(block);
  } //block
  Estimators->stop(estimatorClones);
  for (int ip = 0; ip < NumThreads; ++ip)
    CSMovers[ip]->stopRun2();
#if !defined(REMOVE_TRACEMANAGER)
  Traces->stopRun();
#endif
  //copy back the random states
  for (int ip = 0; ip < NumThreads; ++ip)
    RandomNumberControl::Children[ip] = Rng[ip]->makeClone();
  ///write samples to a file
  bool wrotesamples = DumpConfig;
  if (DumpConfig)
  {
    wrotesamples = MCWalkerConfiguration::dumpEnsemble(wClones, *wOut, myComm->size(), nBlocks);
    if (wrotesamples)
      app_log() << "  samples are written to the config.h5" << std::endl;
  }
  //finalize a qmc section
  return finalize(nBlocks, !wrotesamples);
}

void CSVMC::resetRun()
{
  H1[0]->setPrimary(true);
  IndexType nPsi = Psi1.size();
  for (int ipsi = 1; ipsi < nPsi; ipsi++)
    H1[ipsi]->setPrimary(false);


  makeClones(W, Psi1, H1);
  FairDivideLow(W.getActiveWalkers(), NumThreads, wPerRank);

  if (NumThreads > 1)
    APP_ABORT("OpenMP Parallelization for CSVMC not working at the moment");

  app_log() << "  Initial partition of walkers ";
  copy(wPerRank.begin(), wPerRank.end(), std::ostream_iterator<int>(app_log(), " "));
  app_log() << std::endl;


  if (Movers.empty())
  {
    CSMovers.resize(NumThreads);
    estimatorClones.resize(NumThreads, nullptr);
    traceClones.resize(NumThreads, nullptr);
    Rng.resize(NumThreads);

    // hdf_archive::hdf_archive() is not thread-safe
    for (int ip = 0; ip < NumThreads; ++ip)
      estimatorClones[ip] = new EstimatorManagerBase(*Estimators);

#pragma omp parallel for
    for (int ip = 0; ip < NumThreads; ++ip)
    {
      std::ostringstream os;
      estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
      estimatorClones[ip]->setCollectionMode(false);
#if !defined(REMOVE_TRACEMANAGER)
      traceClones[ip] = Traces->makeClone();
#endif
      Rng[ip] = RandomNumberControl::Children[ip]->makeClone();
      if (qmc_driver_mode[QMC_UPDATE_MODE])
      {
        if (UseDrift == "yes")
        {
          os << "  Using particle-by-particle update with drift " << std::endl;
          CSMovers[ip] = std::make_unique<CSVMCUpdatePbyPWithDriftFast>(*wClones[ip], PsiPoolClones[ip], HPoolClones[ip], *Rng[ip]);
        }
        else
        {
          os << "  Using particle-by-particle update with no drift" << std::endl;
          CSMovers[ip] = std::make_unique<CSVMCUpdatePbyP>(*wClones[ip], PsiPoolClones[ip], HPoolClones[ip], *Rng[ip]);
        }
      }
      else
      {
        if (UseDrift == "yes")
        {
          os << "  Using walker-by-walker update with Drift " << std::endl;
          CSMovers[ip] = std::make_unique<CSVMCUpdateAllWithDrift>(*wClones[ip], PsiPoolClones[ip], HPoolClones[ip], *Rng[ip]);
        }
        else
        {
          os << "  Using walker-by-walker update " << std::endl;
          CSMovers[ip] = std::make_unique<CSVMCUpdateAll>(*wClones[ip], PsiPoolClones[ip], HPoolClones[ip], *Rng[ip]);
        }
      }
      if (ip == 0)
        app_log() << os.str() << std::endl;

      CSMovers[ip]->put(qmcNode);
      CSMovers[ip]->resetRun(branchEngine.get(), estimatorClones[ip], traceClones[ip], DriftModifier);
    }
  }
#if !defined(REMOVE_TRACEMANAGER)
  else
  {
#pragma omp parallel for
    for (int ip = 0; ip < NumThreads; ++ip)
    {
      traceClones[ip]->transfer_state_from(*Traces);
    }
  }
#endif

  app_log() << "  Total Sample Size   =" << nTargetSamples << std::endl;
  app_log() << "  Walker distribution on root = ";
  copy(wPerRank.begin(), wPerRank.end(), std::ostream_iterator<int>(app_log(), " "));
  app_log() << std::endl;
  app_log().flush();

#pragma omp parallel for
  for (int ip = 0; ip < NumThreads; ++ip)
  {
    //int ip=omp_get_thread_num();
    CSMovers[ip]->put(qmcNode);
    CSMovers[ip]->resetRun(branchEngine.get(), estimatorClones[ip], traceClones[ip], DriftModifier);
    if (qmc_driver_mode[QMC_UPDATE_MODE])
      CSMovers[ip]->initCSWalkersForPbyP(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1], nWarmupSteps > 0);
    else
      CSMovers[ip]->initCSWalkers(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1], nWarmupSteps > 0);

    for (int prestep = 0; prestep < nWarmupSteps; ++prestep)
    {
      CSMovers[ip]->advanceWalkers(W.begin() + wPerRank[ip], W.begin() + wPerRank[ip + 1], true);
      if (prestep == nWarmupSteps - 1)
        CSMovers[ip]->updateNorms();
    }
  }

  for (int ip = 0; ip < NumThreads; ++ip)
    wClones[ip]->clearEnsemble();
  if (nSamplesPerThread)
    for (int ip = 0; ip < NumThreads; ++ip)
      wClones[ip]->setNumSamples(nSamplesPerThread);
  nWarmupSteps = 0;
}

} // namespace qmcplusplus

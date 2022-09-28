//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "VMC_CUDA.h"
#include "RandomNumberControl.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/DriftOperators.h"
#include "Utilities/RunTimeManager.h"
#include "Utilities/qmc_common.h"
#ifdef USE_NVTX_API
#include <nvToolsExt.h>
#endif

namespace qmcplusplus
{
/// Constructor.
VMCcuda::VMCcuda(MCWalkerConfiguration& w,
                 TrialWaveFunction& psi,
                 QMCHamiltonian& h,
                 Communicate* comm,
                 bool enable_profiling)
    : QMCDriver(w, psi, h, comm, "VMCcuda", enable_profiling),
      UseDrift("yes"),
      myPeriod4WalkerDump(0)
{
  RootName = "vmc";
  qmc_driver_mode.set(QMC_UPDATE_MODE, 1);
  qmc_driver_mode.set(QMC_WARMUP, 0);
  m_param.add(UseDrift, "useDrift");
  m_param.add(UseDrift, "usedrift");
  m_param.add(nTargetSamples, "targetWalkers");

  H.setRandomGenerator(&Random);
}

bool VMCcuda::checkBounds(std::vector<PosType>& newpos, std::vector<bool>& valid)
{
  for (int iw = 0; iw < newpos.size(); iw++)
  {
    PosType red = W.getLattice().toUnit(newpos[iw]);
    valid[iw]   = W.getLattice().isValid(red);
  }
  return true;
}

// #define DEBUG_DELAYED
// #define SPLIT_SPLINE_DEBUG

void VMCcuda::advanceWalkers()
{
  int nat      = W.getTotalNum();
  int nw       = W.getActiveWalkers();
  int kd       = W.getkblocksize();
  int rcounter = 0;
  std::vector<PosType> delpos(nw);
  std::vector<PosType> newpos(nw);
  std::vector<ValueType> ratios(kd * nw);
  std::vector<RealType> acc_random_nrs(kd * nw);
#ifdef QMC_COMPLEX
  std::vector<RealType> ratios_real(kd * nw);
#endif
  std::vector<GradType> newG(nw * kd);
  std::vector<ValueType> newL(nw * kd);
  std::vector<Walker_t*> accepted(nw);

  for (int isub = 0; isub < nSubSteps; isub++)
  {
#ifdef SPLIT_SPLINE_DEBUG
    if (gpu::rank == 1)
      std::cerr << "sub step: " << isub << "\n";
#endif
    for (int iat = 0; iat < nat; ++iat)
    {
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(delpos, Random);
      for (int iw = 0; iw < nw; ++iw)
      {
        GradType G = W[iw]->G[iat];
        newpos[iw] = W[iw]->R[iat] + m_sqrttau * delpos[iw];
        for (int k = 0; k < kd; k++)
        {
          ratios[iw + k * nw] = 1.0;
#ifdef QMC_COMPLEX
          ratios_real[iw + k * nw] = 1.0;
#endif
        }
      }
      // For each walker, do
      // 1. generate k proposed new positions   [delpos, newpos] -> Done (AT)
      // 2. generate the matrix Rnew (nxk) that contains k proposed new positions -> Done, extended Rnew vector to n vectors (AT)
      W.proposeMove_GPU(newpos, iat);

      // Store acceptance random numbers
      for (int iw = 0; iw < nw; ++iw)
        acc_random_nrs[rcounter++] = Random();

      // check if we are ready to update and if not (k-delay) skip ahead in iat-loop
      if (!W.update_now(iat))
        continue;
      kd           = W.getkupdate();
      int curr_iat = iat - kd + 1;

      Psi.ratio(W, curr_iat, ratios, newG, newL);
      rcounter = 0;
      for (int k = 0; k < kd; ++k)
      {
#ifdef DEBUG_DELAYED
        if (k > 0)
          fprintf(stderr, "*** Delayed Update Path ***\n");
        fprintf(stderr, "iat = %i, k = %i:\n", curr_iat + k, k);
#endif
        accepted.clear();
        for (int iw = 0; iw < nw; ++iw)
          newpos[iw] = W.Rnew_host[iw + k * nw];
        std::vector<bool> acc(nw, true);
        checkBounds(newpos, acc);
#ifdef SPLIT_SPLINE_DEBUG
        if (gpu::rank == 1)
          std::cerr << "iat = " << iat << "\n";
#endif
        if (kDelay)
          Psi.det_lookahead(W, ratios, newG, newL, curr_iat, k, W.getkblocksize(), nw);
#ifdef QMC_COMPLEX
        Psi.convertRatiosFromComplexToReal(ratios, ratios_real);
#endif
#ifdef DEBUG_DELAYED
        for (int iw = 0; iw < nw; ++iw)
          fprintf(stderr, "walker #%i: %f | (%f, %f, %f)\n", iw, ratios[iw + k * nw], newG[iw + k * nw][0],
                  newG[iw + k * nw][1], newG[iw + k * nw][2]);
#endif
        for (int iw = 0; iw < nw; ++iw)
        {
#if defined(DEBUG_DELAYED) || defined(SPLIT_SPLINE_DEBUG)
#ifdef SPLIT_SPLINE_DEBUG
          if (gpu::rank == 1)
#endif
            fprintf(stderr, "Walker #%i (ratio = %f) move ", iw, ratios[iw + k * nw]);
#endif
#ifdef QMC_COMPLEX
          if (acc[iw] && ratios_real[iw + k * nw] * ratios_real[iw + k * nw] > acc_random_nrs[iw + k * nw])
#else
          if (acc[iw] && ratios[iw + k * nw] * ratios[iw + k * nw] > acc_random_nrs[iw + k * nw])
#endif
          {
#ifdef DEBUG_DELAYED
            fprintf(stderr, "accepted.\n");
#endif
            accepted.push_back(W[iw]);
            nAccept++;
            W[iw]->R[curr_iat + k] = W.Rnew_host[iw + k * nw];
            acc[iw]                = true;
          }
          else
          {
#ifdef DEBUG_DELAYED
            fprintf(stderr, "not accepted.\n");
#endif
            acc[iw] = false;
            nReject++;
          }
        }
        W.acceptMove_GPU(acc, k);
        if (accepted.size() || (kDelay > 1))
          Psi.update(&W, accepted, curr_iat, &acc, k);
      }
    }
  }
}

bool VMCcuda::run()
{
  if (UseDrift == "yes")
    return runWithDrift();
#ifdef USE_NVTX_API
  nvtxRangePushA("VMC:run");
#endif
  resetRun();
  IndexType block        = 0;
  IndexType updatePeriod = (qmc_driver_mode[QMC_UPDATE_MODE]) ? Period4CheckProperties : (nBlocks + 1) * nSteps;
  int nat                = W.getTotalNum();
  int nw                 = W.getActiveWalkers();
  std::vector<RealType> LocalEnergy(nw);
  Matrix<ValueType> lapl(nw, nat);
  Matrix<GradType> grad(nw, nat);

  LoopTimer<> vmc_loop;
  RunTimeControl<> runtimeControl(run_time_manager, MaxCPUSecs, myComm->getName(), myComm->rank() == 0);

  // First do warmup steps
  for (int step = 0; step < nWarmupSteps; step++)
  {
#ifdef SPLIT_SPLINE_DEBUG
    if (gpu::rank == 1)
      std::cerr << "Before advanceWalkers(), step " << step << "\n";
#endif
    advanceWalkers();
  }

  // Then accumulate statistics
  do
  {
    vmc_loop.start();
    IndexType step = 0;
    nAccept = nReject = 0;
    Estimators->startBlock(nSteps);
    do
    {
      ++step;
      ++CurrentStep;
      W.resetCollectables();
#ifdef SPLIT_SPLINE_DEBUG
      if (gpu::rank == 1)
        std::cerr << "Before advanceWalkers() loop, step " << step << "\n";
#endif
      advanceWalkers();
      Psi.gradLapl(W, grad, lapl);
      H.evaluate(W, LocalEnergy);
#ifdef SPLIT_SPLINE_DEBUG
      if (gpu::rank == 1)
        for (int ip = 0; ip < nw; ip++)
          fprintf(stderr, "walker #%i energy = %f\n", ip, LocalEnergy[ip]);
#endif
      if (myPeriod4WalkerDump && (CurrentStep % myPeriod4WalkerDump) == 0)
        W.saveEnsemble();
      Estimators->accumulate(W);
    } while (step < nSteps);
    if (nBlocksBetweenRecompute && (1 + block) % nBlocksBetweenRecompute == 0)
      Psi.recompute(W);
    // std::vector<RealType> logPsi(W.WalkerList.size(), 0.0);
    // Psi.evaluateLog(W, logPsi);
    double accept_ratio = (double)nAccept / (double)(nAccept + nReject);
    Estimators->stopBlock(accept_ratio);
    ++block;
    recordBlock(block);
    vmc_loop.stop();

    bool stop_requested = false;
    // Rank 0 decides whether the time limit was reached
    if (!myComm->rank())
      stop_requested = runtimeControl.checkStop(vmc_loop);
    myComm->bcast(stop_requested);
    if (stop_requested)
    {
      if (!myComm->rank())
        app_log() << runtimeControl.generateStopMessage("VMC_CUDA", block);
      run_time_manager.markStop();
      break;
    }

  } while (block < nBlocks);
  //Mover->stopRun();
  //finalize a qmc section
  if (!myComm->rank())
  {
    std::cerr << "At the end of VMC" << std::endl;
    gpu::cuda_memory_manager.report();
  }
#ifdef SPLIT_SPLINE_DEBUG
  _exit(0);
#endif
#ifdef USE_NVTX_API
  nvtxRangePop();
#endif
  return finalize(block);
}

void VMCcuda::advanceWalkersWithDrift()
{
  bool update_now = false;
  int nat         = W.getTotalNum();
  int nw          = W.getActiveWalkers();
  std::vector<PosType> delpos(nw);
  std::vector<PosType> newpos(nw);
  std::vector<ValueType> ratios(nw);
#ifdef QMC_COMPLEX
  std::vector<RealType> ratios_real(nw);
#endif
  std::vector<GradType> oldG(nw), newG(nw);
  std::vector<ValueType> oldL(nw), newL(nw);
  std::vector<Walker_t*> accepted(nw);
  std::vector<bool> acc(nw, true);

  for (int isub = 0; isub < nSubSteps; isub++)
  {
    int k = 0;
    for (int iat = 0; iat < nat; iat++)
    {
      Psi.calcGradient(W, iat, k, oldG);
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(delpos, Random);
      Psi.addGradient(W, iat, oldG);
      for (int iw = 0; iw < nw; iw++)
      {
        PosType dr;
        delpos[iw] *= m_sqrttau;
        DriftModifier->getDrift(m_tauovermass, oldG[iw], dr);
        newpos[iw] = W[iw]->R[iat] + delpos[iw] + dr;
        ratios[iw] = 1.0;
        acc[iw]    = true;
#ifdef QMC_COMPLEX
        ratios_real[iw] = 1.0;
#endif
      }
      W.proposeMove_GPU(newpos, iat);
      update_now = W.update_now(iat);
      Psi.calcRatio(W, iat, ratios, newG, newL);
      accepted.clear();
      checkBounds(newpos, acc);
      if (kDelay)
        Psi.det_lookahead(W, ratios, newG, newL, iat, k, W.getkblocksize(), nw);
      std::vector<RealType> logGf_v(nw);
      std::vector<RealType> rand_v(nw);
      for (int iw = 0; iw < nw; ++iw)
      {
        logGf_v[iw] = -m_oneover2tau * dot(delpos[iw], delpos[iw]);
        rand_v[iw]  = Random();
      }
      Psi.addRatio(W, iat, k, ratios, newG, newL);
#ifdef QMC_COMPLEX
      Psi.convertRatiosFromComplexToReal(ratios, ratios_real);
#endif
      for (int iw = 0; iw < nw; ++iw)
      {
        PosType drNew;
        DriftModifier->getDrift(m_tauovermass, newG[iw], drNew);
        drNew += newpos[iw] - W[iw]->R[iat];
        // if (dot(drNew, drNew) > 25.0)
        //   std::cerr << "Large drift encountered!  Drift = " << drNew << std::endl;
        RealType logGb = -m_oneover2tau * dot(drNew, drNew);
        RealType x     = logGb - logGf_v[iw];
#ifdef QMC_COMPLEX
        RealType prob = ratios_real[iw] * ratios_real[iw] * std::exp(x);
#else
        RealType prob = ratios[iw] * ratios[iw] * std::exp(x);
#endif
#ifdef DEBUG_DELAYED
        fprintf(stderr, " Walker #%i (ratio = %f) move ", iw, ratios[iw]);
#endif
        if (acc[iw] && rand_v[iw] < prob)
        {
#ifdef DEBUG_DELAYED
          fprintf(stderr, "accepted.\n");
#endif
          accepted.push_back(W[iw]);
          nAccept++;
          W[iw]->R[iat] = newpos[iw];
          acc[iw]       = true;
        }
        else
        {
#ifdef DEBUG_DELAYED
          fprintf(stderr, "not accepted.\n");
#endif
          acc[iw] = false;
          nReject++;
        }
      }
      W.acceptMove_GPU(acc, k);
      if (accepted.size() || (kDelay > 1))
        Psi.update(&W, accepted, iat, &acc, k);
      if (kDelay > 1)
        k++;
      if (update_now)
        k = 0;
    }
  }
}

bool VMCcuda::runWithDrift()
{
#ifdef USE_NVTX_API
  nvtxRangePushA("VMC:runWithDrift");
#endif
  resetRun();
  IndexType block = 0;
  int nat         = W.getTotalNum();
  int nw          = W.getActiveWalkers();
  std::vector<RealType> LocalEnergy(nw);
  Matrix<ValueType> lapl(nw, nat);
  Matrix<GradType> grad(nw, nat);

  // First, do warmup steps
  for (int step = 0; step < nWarmupSteps; step++)
    advanceWalkersWithDrift();

  // Now do data collection steps
  do
  {
    IndexType step = 0;
    nAccept = nReject = 0;
    Estimators->startBlock(nSteps);
    do
    {
      step++;
      CurrentStep++;
      W.resetCollectables();
      advanceWalkersWithDrift();
      Psi.gradLapl(W, grad, lapl);
      H.evaluate(W, LocalEnergy);
      if (myPeriod4WalkerDump && (CurrentStep % myPeriod4WalkerDump) == 0)
        W.saveEnsemble();
      Estimators->accumulate(W);
    } while (step < nSteps);
    if (nBlocksBetweenRecompute && (1 + block) % nBlocksBetweenRecompute == 0)
      Psi.recompute(W);
    //      Psi.recompute(W);
    double accept_ratio = (double)nAccept / (double)(nAccept + nReject);
    Estimators->stopBlock(accept_ratio);
    ++block;
    recordBlock(block);
  } while (block < nBlocks);
  //finalize a qmc section
  if (!myComm->rank())
  {
    std::cerr << "At the end of VMC with drift" << std::endl;
    gpu::cuda_memory_manager.report();
  }
#ifdef USE_NVTX_API
  nvtxRangePop();
#endif
  return finalize(block);
}


void VMCcuda::resetRun()
{
  SpeciesSet tspecies(W.getSpeciesSet());
  int massind              = tspecies.addAttribute("mass");
  RealType mass            = tspecies(massind, 0);
  RealType oneovermass     = 1.0 / mass;
  RealType oneoversqrtmass = std::sqrt(oneovermass);
  m_oneover2tau            = 0.5 * mass / Tau;
  m_sqrttau                = std::sqrt(Tau / mass);
  m_tauovermass            = Tau / mass;
  if (!myComm->rank())
  {
    std::cerr << "Before allocating GPU buffer" << std::endl;
    gpu::cuda_memory_manager.report();
  }
  // Compute the size of data needed for each walker on the GPU card
  PointerPool<Walker_t::cuda_Buffer_t> pool;
  app_log() << "Starting VMCcuda::resetRun() " << std::endl;
  Psi.reserve(pool, false, W.getkblocksize());
  app_log() << "Each walker requires " << pool.getTotalSize() * sizeof(CTS::ValueType) << " bytes in GPU memory.\n";
  // Now allocate memory on the GPU card for each walker
  // for (int iw=0; iw<W.WalkerList.size(); iw++) {
  //   Walker_t &walker = *(W.WalkerList[iw]);
  //   walker.resizeCuda(pool.getTotalSize());
  //   // pool.allocate(walker.cuda_DataSet);
  // }
  W.allocateGPU(pool.getTotalSize());
  if (!myComm->rank())
  {
    std::cerr << "After allocating GPU buffer" << std::endl;
    gpu::cuda_memory_manager.report();
  }
  app_log() << "Successfully allocated walkers.\n";
  W.copyWalkersToGPU();
  W.updateLists_GPU();
  std::vector<RealType> logPsi(W.WalkerList.size(), 0.0);
  //Psi.evaluateLog(W, logPsi);
  Psi.recompute(W, true);
  Estimators->start(nBlocks, true);
  // Compute sample dumping frequency
  if (nTargetSamples > (myComm->size() * W.WalkerList.size()))
  {
    int nnodes           = myComm->size();
    int nw               = W.WalkerList.size();
    int samples_per_node = (nTargetSamples + nnodes - 1) / nnodes;
    int dumps_per_node   = (samples_per_node + nw - 1) / nw;
    myPeriod4WalkerDump  = Period4WalkerDump;
    app_log() << "  Dumping walker ensemble every " << myPeriod4WalkerDump << " steps.\n";
  }
  else
  {
    myPeriod4WalkerDump = nBlocks * nSteps;
    nTargetSamples      = myComm->size() * W.WalkerList.size();
  }
  W.clearEnsemble();
  int samples_this_node = nTargetSamples / myComm->size();
  if (nTargetSamples % myComm->size() > myComm->rank())
    samples_this_node += 1;
  app_log() << "  Node zero will generate " << samples_this_node << " samples.\n";
  W.setNumSamples(samples_this_node);
}

bool VMCcuda::put(xmlNodePtr q)
{
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
    int Nthreads = 1;
    int Nprocs   = myComm->size();
    //target samples set by samples or samplesperthread/dmcwalkersperthread
    nTargetPopulation = std::max(nTargetPopulation, nSamplesPerThread * Nprocs * Nthreads);
    nTargetSamples    = static_cast<int>(std::ceil(nTargetPopulation));
    if (nTargetSamples)
    {
      int nwtot      = nw * Nprocs; //total number of walkers used by this qmcsection
      nTargetSamples = std::max(nwtot, nTargetSamples);
      nTargetSamples = ((nTargetSamples + nwtot - 1) / nwtot) *
          nwtot; // nTargetSamples are always multiples of total number of walkers
      nSamplesPerThread = nTargetSamples / Nprocs / Nthreads;
      int ns_target     = nTargetSamples * nStepsBetweenSamples; //total samples to generate
      int ns_per_step   = Nprocs * nw;                           //total samples per step
      nSteps            = std::max(nSteps, (ns_target / ns_per_step + nBlocks - 1) / nBlocks);
      Period4WalkerDump = nStepsBetweenSamples = ns_per_step * nSteps * nBlocks / nTargetSamples;
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

  if (DumpConfig)
    app_log() << "  DumpConfig==true Configurations are dumped to config.h5 with a period of " << Period4CheckPoint
              << " blocks" << std::endl;
  else
    app_log() << "  DumpConfig==false Nothing (configurations, state) will be saved." << std::endl;
  if (Period4WalkerDump > 0)
    app_log() << "  Walker Samples are dumped every " << Period4WalkerDump << " steps." << std::endl;
  app_log() << "</vmc>" << std::endl;
  app_log().flush();
  //nothing to add
  return true;
}

} // namespace qmcplusplus

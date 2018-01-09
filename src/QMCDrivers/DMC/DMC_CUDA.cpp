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
    
    


#include "QMCDrivers/DMC/DMC_CUDA.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/QMCUpdateBase.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "Utilities/RunTimeManager.h"
#include "Message/CommOperators.h"
#include "type_traits/scalar_traits.h"


namespace qmcplusplus
{

/// Constructor.
DMCcuda::DMCcuda(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                 QMCHamiltonian& h,WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool), myWarmupSteps(0), Mover(0),
  NLop(w.getTotalNum()),
  ResizeTimer("DMCcuda::resize"),
  DriftDiffuseTimer("DMCcuda::Drift_Diffuse"),
  BranchTimer("DMCcuda::Branch"),
  HTimer("DMCcuda::Hamiltonian")
{
  RootName = "dmc";
  QMCType ="DMCcuda";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  //m_param.add(myWarmupSteps,"warmupSteps","int");
  //m_param.add(nTargetSamples,"targetWalkers","int");
  m_param.add(ScaleWeight, "scaleweight", "string");
  TimerManager.addTimer (&ResizeTimer);
  TimerManager.addTimer (&DriftDiffuseTimer);
  TimerManager.addTimer (&BranchTimer);
  TimerManager.addTimer (&HTimer);
}

bool DMCcuda::checkBounds (const PosType &newpos)
{
  if (!W.UseBoundBox)
    return true;
  PosType red = W.Lattice.toUnit(newpos);
  return W.Lattice.isValid(red);
}

void DMCcuda::checkBounds (std::vector<PosType> &newpos,
                           std::vector<bool> &valid)
{
  for (int iw=0; iw<newpos.size(); iw++)
  {
    PosType red = W.Lattice.toUnit(newpos[iw]);
    valid[iw] = W.Lattice.isValid(red);
  }
}


bool DMCcuda::run()
{
  bool scaleweight = ScaleWeight == "yes";
  if (scaleweight)
    app_log() << "  Scaling weight per Umrigar/Nightingale.\n";
  resetRun();
  Mover->MaxAge = 1;
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  int nat = W.getTotalNum();
  int nw  = W.getActiveWalkers();
  std::vector<RealType>  LocalEnergy(nw), LocalEnergyOld(nw);
  std::vector<PosType>   delpos(nw);
  std::vector<PosType>   newpos(nw);
  std::vector<ValueType> ratios(nw), rplus(nw), rminus(nw);
  std::vector<RealType>  R2prop(nw), R2acc(nw);
#ifdef QMC_COMPLEX
  std::vector<RealType>  ratios_real(nw);
#endif
  std::vector<GradType>  oldG(nw), newG(nw);
  std::vector<ValueType> oldL(nw), newL(nw);
  std::vector<Walker_t*> accepted(nw);
  Matrix<ValueType> lapl(nw, nat);
  Matrix<GradType>  grad(nw, nat);
  std::vector<std::vector<NonLocalData> > Txy(nw);
  for (int iw=0; iw<nw; iw++)
    W[iw]->Weight = 1.0;

  LoopTimer dmc_loop;
  RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);
  bool enough_time_for_next_iteration = true;
  do
  {
    dmc_loop.start();
    IndexType step = 0;
    nAccept = nReject = 0;
    Estimators->startBlock(nSteps);
    do
    {
      step++;
      CurrentStep++;
      nw = W.getActiveWalkers();
      ResizeTimer.start();
      LocalEnergy.resize(nw);
      delpos.resize(nw);
      newpos.resize(nw);
      ratios.resize(nw);
#ifdef QMC_COMPLEX
      ratios_real.resize(nw);
#endif
      rplus.resize(nw);
      rminus.resize(nw);
      oldG.resize(nw);
      newG.resize(nw);
      oldL.resize(nw);
      newL.resize(nw);
      accepted.resize(nw);
      lapl.resize(nw, nat);
      grad.resize(nw, nat);
      R2prop.resize(nw,0.0);
      R2acc.resize(nw,0.0);
      W.updateLists_GPU();
      ResizeTimer.stop();
      if (UseTMove)
      {
        Txy.resize(nw);
        for (int iw=0; iw<nw; iw++)
        {
          Txy[iw].clear();
          Txy[iw].push_back(NonLocalData(-1, 1.0, PosType()));
        }
      }
      for (int iw=0; iw<nw; iw++)
        W[iw]->Age++;
      DriftDiffuseTimer.start();
      for(int iat=0; iat<nat; iat++)
      {
        Psi.calcGradient (W, iat, oldG);
        //create a 3N-Dimensional Gaussian with variance=1
        makeGaussRandomWithEngine(delpos,Random);
        Psi.addGradient(W, iat, oldG);
        for(int iw=0; iw<nw; iw++)
        {
          delpos[iw] *= m_sqrttau;
          PosType dr;
          getScaledDrift(m_tauovermass,oldG[iw],dr);
          newpos[iw] = W[iw]->R[iat] + delpos[iw] + dr;
          ratios[iw] = 1.0;
#ifdef QMC_COMPLEX
          ratios_real[iw] = 1.0;
#endif
          R2prop[iw] += dot(delpos[iw], delpos[iw]);
        }
        W.proposeMove_GPU(newpos, iat);
        Psi.calcRatio(W,iat,ratios,newG, newL);
        accepted.clear();
        std::vector<bool> acc(nw, true);
        if (W.UseBoundBox)
          checkBounds (newpos, acc);
        std::vector<RealType> logGf_v(nw);
        std::vector<RealType> rand_v(nw);
        for(int iw=0; iw<nw; ++iw)
        {
          logGf_v[iw] = -m_oneover2tau * dot(delpos[iw], delpos[iw]);
          rand_v[iw] = Random();
        }
        Psi.addRatio(W, iat, ratios, newG, newL);
#ifdef QMC_COMPLEX
        Psi.convertRatiosFromComplexToReal(ratios, ratios_real);
#endif
        for(int iw=0; iw<nw; ++iw)
        {
          PosType drNew;
          getScaledDrift(m_tauovermass,newG[iw],drNew);
          drNew += newpos[iw] - W[iw]->R[iat];
          RealType logGb = -m_oneover2tau * dot(drNew, drNew);
          RealType x = logGb - logGf_v[iw];

#ifdef QMC_COMPLEX
          RealType prob = ratios_real[iw]*ratios_real[iw]*std::exp(x);
          if(acc[iw] && rand_v[iw] < prob && ratios_real[iw] > 0.0)
#else
          RealType prob = ratios[iw]*ratios[iw]*std::exp(x);
          if(acc[iw] && rand_v[iw] < prob && ratios[iw] > 0.0)
#endif
          {
            accepted.push_back(W[iw]);
            nAccept++;
            W[iw]->R[iat] = newpos[iw];
            W[iw]->Age = 0;
            acc[iw] = true;
            R2acc[iw] += dot(delpos[iw], delpos[iw]);
          }
          else
          {
            acc[iw] = false;
            nReject++;
          }
        }
        W.acceptMove_GPU(acc);
        if (accepted.size())
          Psi.update(accepted,iat);
      }
      DriftDiffuseTimer.stop();
      Psi.gradLapl(W, grad, lapl);
      HTimer.start();
      if (UseTMove)
        H.evaluate (W, LocalEnergy, Txy);
      else
        H.evaluate (W, LocalEnergy);
      HTimer.stop();
      if (CurrentStep == 1)
        LocalEnergyOld = LocalEnergy;
      if (UseTMove==TMOVE_V0)
      {
        // Now, attempt nonlocal move
        accepted.clear();
        std::vector<int> iatList;
        std::vector<PosType> accPos;
        for (int iw=0; iw<nw; iw++)
        {
          int ibar = NLop.selectMove(Random(), Txy[iw]);
          if (ibar)
          {
            int iat = Txy[iw][ibar].PID;
            PosType newpos(W[iw]->R[iat] + Txy[iw][ibar].Delta);
            if (checkBounds(newpos))
            {
              accepted.push_back(W[iw]);
              iatList.push_back(iat);
              accPos.push_back(newpos);
            }
          }
        }
        if (accepted.size())
        {
          Psi.ratio(accepted,iatList, accPos, ratios, newG, newL);
          Psi.update(accepted,iatList);
          for (int i=0; i<accepted.size(); i++)
            accepted[i]->R[iatList[i]] = accPos[i];
          W.NLMove_GPU (accepted, accPos, iatList);
        }
      }
      else if(UseTMove==TMOVE_V1)
      {
        APP_ABORT("Tmove v1 has not been implemented on GPU.\n  please contact the developers if you need this feature");
      }
      // Now branch
      BranchTimer.start();
      for (int iw=0; iw<nw; iw++)
      {
        RealType v2=0.0, v2bar=0.0;
        for(int iat=0; iat<nat; iat++)
        {
          PosType wG_scaled;
          getScaledDrift(m_tauovermass,W.G[iat],wG_scaled);
          v2bar += dot(wG_scaled,wG_scaled);
#ifdef QMC_COMPLEX
          PosType wG_real;
          convert(W.G[iat],wG_real);
          v2 += dot(wG_real,wG_real);
#else
          // should be removed when things work fine
          v2 += dot(W.G[iat],W.G[iat]);
#endif
        }
        RealType scNew = std::sqrt(v2bar/(v2*m_tauovermass*m_tauovermass));
        RealType scOld = (CurrentStep == 1) ? scNew : W[iw]->getPropertyBase()[DRIFTSCALE];
        W[iw]->getPropertyBase()[DRIFTSCALE] = scNew;
        // fprintf (stderr, "iw = %d  scNew = %1.8f  scOld = %1.8f\n", iw, scNew, scOld);
        RealType tauRatio = R2acc[iw] / R2prop[iw];
        //allow large time steps during warmup
        //if (tauRatio < 0.5)
        //  std::cerr << "  tauRatio = " << tauRatio << std::endl;
        RealType taueff = m_tauovermass * tauRatio;
        if (scaleweight)
          W[iw]->Weight *= branchEngine->branchWeightTau
                           (LocalEnergy[iw], LocalEnergyOld[iw], scNew, scOld, taueff);
        else
          W[iw]->Weight *= branchEngine->branchWeight
                           (LocalEnergy[iw], LocalEnergyOld[iw]);
        W[iw]->getPropertyBase()[R2ACCEPTED] = R2acc[iw];
        W[iw]->getPropertyBase()[R2PROPOSED] = R2prop[iw];
      }
      Mover->setMultiplicity(W.begin(), W.end());
      branchEngine->branch(CurrentStep,W);
      nw = W.getActiveWalkers();
      LocalEnergyOld.resize(nw);
      for (int iw=0; iw<nw; iw++)
        LocalEnergyOld[iw] = W[iw]->getPropertyBase()[LOCALENERGY];
      BranchTimer.stop();
    }
    while(step<nSteps);
    if ( nBlocksBetweenRecompute && (1+block)%nBlocksBetweenRecompute == 0 ) Psi.recompute(W, true);
    double accept_ratio = (double)nAccept/(double)(nAccept+nReject);
    Estimators->stopBlock(accept_ratio);
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    ++block;
    recordBlock(block);
    dmc_loop.stop();
    enough_time_for_next_iteration = runtimeControl.enough_time_for_next_iteration(dmc_loop);
    // Rank 0 decides whether the time limit was reached
    myComm->bcast(enough_time_for_next_iteration);
    if (!enough_time_for_next_iteration)
    {
      app_log() << runtimeControl.time_limit_message("DMC", block);
    }
  }
  while(block<nBlocks && enough_time_for_next_iteration);
  //finalize a qmc section
  return finalize(block);
}


void DMCcuda::resetUpdateEngine()
{
  if(Mover==0) //disable switching update modes for DMC in a run
  {
    //load walkers if they were saved
    W.loadEnsemble();
    branchEngine->initWalkerController(W,false,false);
    Mover = new DMCUpdatePbyPWithRejectionFast(W,Psi,H,Random);
    Mover->resetRun(branchEngine,Estimators);
    //Mover->initWalkersForPbyP(W.begin(),W.end());
  }
  else
  {
    int nw_multi=branchEngine->resetRun(qmcNode);
    if(nw_multi>1)
    {
      app_log() << " Current population " << W.getActiveWalkers() << " " << W.getGlobalNumWalkers()  << std::endl;
      app_log() << " The target population has changed. Multiply walkers by " << nw_multi << std::endl;
      W.createWalkers((nw_multi-1)*W.getActiveWalkers());
      setWalkerOffsets();
      app_log() << " New population " << W.getActiveWalkers() << " " << W.getGlobalNumWalkers()  << std::endl;
    }
  }

  branchEngine->checkParameters(W);

  //    Mover->updateWalkers(W.begin(),W.end());
  Mover->MaxAge=1;
  app_log() << "  Steps per block = " << nSteps << std::endl;
  app_log() << "  Number of blocks = " << nBlocks << std::endl;
}



void DMCcuda::resetRun()
{
  resetUpdateEngine();
  SpeciesSet tspecies(W.getSpeciesSet());
  int massind=tspecies.addAttribute("mass");
  RealType mass = tspecies(massind,0);
  RealType oneovermass = 1.0/mass;
  RealType oneoversqrtmass = std::sqrt(oneovermass);
  m_oneover2tau = 0.5/Tau;
  m_sqrttau = std::sqrt(Tau/mass);
  m_tauovermass = Tau/mass;
  if (!myComm->rank())
    gpu::cuda_memory_manager.report();
  // Compute the size of data needed for each walker on the GPU card
  PointerPool<Walker_t::cuda_Buffer_t > pool;
  Psi.reserve (pool);
  app_log() << "Each walker requires "
            << pool.getTotalSize() * sizeof(CudaValueType)
            << " bytes in GPU memory.\n";
  app_log() << "Preparing to allocate " << W.WalkerList.size()
            << " walkers.\n";
  // Now allocate memory on the GPU card for each walker
  int cudaSize = pool.getTotalSize();
  // for (int iw=0; iw<W.WalkerList.size(); iw++) {
  //   Walker_t &walker = *(W.WalkerList[iw]);
  //   walker.resizeCuda(cudaSize);
  //   //pool.allocate(walker.cuda_DataSet);
  // }
  W.allocateGPU(pool.getTotalSize());
  app_log() << "Successfully allocated walkers.\n";
  W.copyWalkersToGPU();
  W.updateLists_GPU();
  std::vector<RealType> logPsi(W.WalkerList.size(), 0.0);
  //Psi.evaluateLog(W, logPsi);
  Psi.recompute(W, true);
  Estimators->start(nBlocks, true);
}

bool
DMCcuda::put(xmlNodePtr q)
{
  //nothing to add
  UseTMove = NLop.put(q);

  BranchInterval=-1;
  ParameterSet p;
  p.add(BranchInterval,"branchInterval","string");
  p.add(BranchInterval,"branchinterval","string");
  p.add(BranchInterval,"substeps","int");
  p.add(BranchInterval,"subSteps","int");
  p.add(BranchInterval,"sub_steps","int");
  p.put(q);
  return true;
}
}


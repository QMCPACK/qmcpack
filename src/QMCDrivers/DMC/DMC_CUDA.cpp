//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Ken Esler
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMC/DMC_CUDA.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/QMCUpdateBase.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"


namespace qmcplusplus
{

/// Constructor.
DMCcuda::DMCcuda(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                 QMCHamiltonian& h,WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool), myWarmupSteps(0), Mover(0),
  ResizeTimer("DMCcuda::resize"),
  DriftDiffuseTimer("DMCcuda::Drift/Diffuse"),
  BranchTimer("DMCcuda::Branch"),
  HTimer("DMCcuda::Hamiltonian")
{
  RootName = "dmc";
  QMCType ="DMCcuda";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  //m_param.add(myWarmupSteps,"warmupSteps","int");
  //m_param.add(nTargetSamples,"targetWalkers","int");
  m_param.add(NonLocalMove,"nonlocalmove","string");
  m_param.add(NonLocalMove,"nonlocalmoves","string");
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

void DMCcuda::checkBounds (vector<PosType> &newpos,
                           vector<bool> &valid)
{
  for (int iw=0; iw<newpos.size(); iw++)
  {
    PosType red = W.Lattice.toUnit(newpos[iw]);
    valid[iw] = W.Lattice.isValid(red);
  }
}


bool DMCcuda::run()
{
  bool NLmove = NonLocalMove == "yes";
  bool scaleweight = ScaleWeight == "yes";
  if (NLmove)
    app_log() << "  Using Casula nonlocal moves in DMCcuda.\n";
  if (scaleweight)
    app_log() << "  Scaling weight per Umrigar/Nightengale.\n";
  resetRun();
  Mover->MaxAge = 1;
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  int nat = W.getTotalNum();
  int nw  = W.getActiveWalkers();
  vector<RealType>  LocalEnergy(nw), LocalEnergyOld(nw),
         oldScale(nw), newScale(nw);
  vector<PosType>   delpos(nw);
  vector<PosType>   dr(nw);
  vector<PosType>   newpos(nw);
  vector<ValueType> ratios(nw), rplus(nw), rminus(nw), R2prop(nw), R2acc(nw);
  vector<PosType>  oldG(nw), newG(nw);
  vector<ValueType> oldL(nw), newL(nw);
  vector<Walker_t*> accepted(nw);
  Matrix<ValueType> lapl(nw, nat);
  Matrix<GradType>  grad(nw, nat);
  vector<ValueType> V2(nw), V2bar(nw);
  vector<vector<NonLocalData> > Txy(nw);
  for (int iw=0; iw<nw; iw++)
    W[iw]->Weight = 1.0;
  do
  {
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
      oldScale.resize(nw);
      newScale.resize(nw);
      delpos.resize(nw);
      dr.resize(nw);
      newpos.resize(nw);
      ratios.resize(nw);
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
      V2.resize(nw,0.0);
      V2bar.resize(nw,0.0);
      W.updateLists_GPU();
      ResizeTimer.stop();
      if (NLmove)
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
          oldScale[iw] = getDriftScale(m_tauovermass,oldG[iw]);
          dr[iw] = delpos[iw] + (oldScale[iw]*oldG[iw]);
          newpos[iw]=W[iw]->R[iat] + dr[iw];
          ratios[iw] = 1.0;
          R2prop[iw] += dot(delpos[iw], delpos[iw]);
        }
        W.proposeMove_GPU(newpos, iat);
        Psi.calcRatio(W,iat,ratios,newG, newL);
        accepted.clear();
        vector<bool> acc(nw, true);
        if (W.UseBoundBox)
          checkBounds (newpos, acc);
        std::vector<RealType> logGf_v(nw);
        std::vector<RealType> rand_v(nw);
        for(int iw=0; iw<nw; ++iw)
        {
          PosType drOld =
            newpos[iw] - (W[iw]->R[iat] + oldScale[iw]*oldG[iw]);
          logGf_v[iw] = -m_oneover2tau * dot(drOld, drOld);
          rand_v[iw] = Random();
        }
        Psi.addRatio(W,iat,ratios, newG, newL);
        for(int iw=0; iw<nw; ++iw)
        {
          newScale[iw]   = getDriftScale(m_tauovermass,newG[iw]);
          PosType drNew  =
            (newpos[iw] + newScale[iw]*newG[iw]) - W[iw]->R[iat];
          RealType logGb =  -m_oneover2tau * dot(drNew, drNew);
          RealType x = logGb - logGf_v[iw];
          RealType prob = ratios[iw]*ratios[iw]*std::exp(x);
          if(acc[iw] && rand_v[iw] < prob && ratios[iw] > 0.0)
          {
            accepted.push_back(W[iw]);
            nAccept++;
            W[iw]->R[iat] = newpos[iw];
            W[iw]->Age = 0;
            acc[iw] = true;
            R2acc[iw] += dot(delpos[iw], delpos[iw]);
            V2[iw]    += m_tauovermass * m_tauovermass * dot(newG[iw],newG[iw]);
            V2bar[iw] +=  newScale[iw] *  newScale[iw] * dot(newG[iw],newG[iw]);
          }
          else
          {
            acc[iw] = false;
            nReject++;
            V2[iw]    += m_tauovermass * m_tauovermass * dot(oldG[iw],oldG[iw]);
            V2bar[iw] +=  oldScale[iw] *  oldScale[iw] * dot(oldG[iw],oldG[iw]);
          }
        }
        W.acceptMove_GPU(acc);
        if (accepted.size())
          Psi.update(accepted,iat);
      }
      DriftDiffuseTimer.stop();
      //	Psi.recompute(W, false);
      Psi.gradLapl(W, grad, lapl);
      HTimer.start();
      if (NLmove)	  H.evaluate (W, LocalEnergy, Txy);
      else    	  H.evaluate (W, LocalEnergy);
      HTimer.stop();
// 	for (int iw=0; iw<nw; iw++) {
// 	  branchEngine->clampEnergy(LocalEnergy[iw]);
// 	  W[iw]->getPropertyBase()[LOCALENERGY] = LocalEnergy[iw];
// 	}
      if (CurrentStep == 1)
        LocalEnergyOld = LocalEnergy;
      if (NLmove)
      {
        // Now, attempt nonlocal move
        accepted.clear();
        vector<int> iatList;
        vector<PosType> accPos;
        for (int iw=0; iw<nw; iw++)
        {
          /// HACK HACK HACK
// 	    if (LocalEnergy[iw] < -2300.0) {
// 	      cerr << "Walker " << iw << " has energy "
// 		   << LocalEnergy[iw] << endl;;
// 	      double maxWeight = 0.0;
// 	      int elMax = -1;
// 	      PosType posMax;
// 	      for (int j=1; j<Txy[iw].size(); j++)
// 		if (std::fabs(Txy[iw][j].Weight) > std::fabs(maxWeight)) {
// 		  maxWeight = Txy[iw][j].Weight;
// 		  elMax = Txy[iw][j].PID;
// 		  posMax = W[iw]->R[elMax] + Txy[iw][j].Delta;
// 		}
// 	      cerr << "Maximum weight is " << maxWeight << " for electron "
// 		   << elMax << " at position " << posMax << endl;
// 	      PosType unit = W.Lattice.toUnit(posMax);
// 	      unit[0] -= round(unit[0]);
// 	      unit[1] -= round(unit[1]);
// 	      unit[2] -= round(unit[2]);
// 	      cerr << "Reduced position = " << unit << endl;
// 	    }
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
          // HACK HACK HACK
          // Recompute the kinetic energy
          // Psi.gradLapl(W, grad, lapl);
          // H.evaluate (W, LocalEnergy);
          //W.copyWalkersToGPU();
        }
      }
      // Now branch
      BranchTimer.start();
      for (int iw=0; iw<nw; iw++)
      {
        RealType v2=0.0, v2bar=0.0;
        for(int iat=0; iat<nat; iat++)
        {
          v2 += dot(W.G[iat],W.G[iat]);
          RealType newscale = getDriftScale(m_tauovermass,newG[iw]);
          v2 += m_tauovermass * m_tauovermass * dot(newG[iw],newG[iw]);
          v2bar +=  newscale * newscale * dot(newG[iw],newG[iw]);
        }
        //RealType scNew = std::sqrt(V2bar[iw] / V2[iw]);
        RealType scNew = std::sqrt(v2bar/v2);
        RealType scOld = (CurrentStep == 1) ? scNew : W[iw]->getPropertyBase()[DRIFTSCALE];
        W[iw]->getPropertyBase()[DRIFTSCALE] = scNew;
        // fprintf (stderr, "iw = %d  scNew = %1.8f  scOld = %1.8f\n", iw, scNew, scOld);
        RealType tauRatio = R2acc[iw] / R2prop[iw];
        //allow large time steps during warmup
        //if (tauRatio < 0.5)
        //  cerr << "  tauRatio = " << tauRatio << endl;
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
    Psi.recompute(W, true);
    double accept_ratio = (double)nAccept/(double)(nAccept+nReject);
    Estimators->stopBlock(accept_ratio);
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    ++block;
    recordBlock(block);
  }
  while(block<nBlocks);
  //finalize a qmc section
  return finalize(block);
}



bool DMCcuda::runWithNonlocal()
{
  resetRun();
  Mover->MaxAge = 1;
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  int nat = W.getTotalNum();
  int nw  = W.getActiveWalkers();
  vector<RealType>  LocalEnergy(nw), LocalEnergyOld(nw),
         oldScale(nw), newScale(nw);
  vector<PosType>   delpos(nw);
  vector<PosType>   dr(nw);
  vector<PosType>   newpos(nw);
  vector<ValueType> ratios(nw), rplus(nw), rminus(nw), R2prop(nw), R2acc(nw);
  vector<PosType>  oldG(nw), newG(nw);
  vector<ValueType> oldL(nw), newL(nw);
  vector<Walker_t*> accepted(nw);
  Matrix<ValueType> lapl(nw, nat);
  Matrix<GradType>  grad(nw, nat);
  vector<vector<NonLocalData> > Txy(nw);
  for (int iw=0; iw<nw; iw++)
    W[iw]->Weight = 1.0;
  do
  {
    IndexType step = 0;
    nAccept = nReject = 0;
    Estimators->startBlock(nSteps);
    do
    {
      step++;
      CurrentStep++;
      nw = W.getActiveWalkers();
      LocalEnergy.resize(nw);
      oldScale.resize(nw);
      newScale.resize(nw);
      delpos.resize(nw);
      dr.resize(nw);
      newpos.resize(nw);
      ratios.resize(nw);
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
      Txy.resize(nw);
      for (int iw=0; iw<nw; iw++)
      {
        Txy[iw].clear();
        Txy[iw].push_back(NonLocalData(-1, 1.0, PosType()));
        W[iw]->Age++;
      }
      for(int iat=0; iat<nat; iat++)
      {
        Psi.calcGradient (W, iat, oldG);
        //create a 3N-Dimensional Gaussian with variance=1
        makeGaussRandomWithEngine(delpos,Random);
        Psi.addGradient(W, iat, oldG);
        for(int iw=0; iw<nw; iw++)
        {
          delpos[iw] *= m_sqrttau;
          oldScale[iw] = getDriftScale(m_tauovermass,oldG[iw]);
          dr[iw] = delpos[iw] + (oldScale[iw]*oldG[iw]);
          newpos[iw]=W[iw]->R[iat] + dr[iw];
          ratios[iw] = 1.0;
          R2prop[iw] += dot(delpos[iw], delpos[iw]);
        }
        W.proposeMove_GPU(newpos, iat);
        Psi.calcRatio(W,iat,ratios,newG, newL);
        accepted.clear();
        vector<bool> acc(nw, false);
        std::vector<RealType> logGf_v(nw);
        std::vector<RealType> rand_v(nw);
        for(int iw=0; iw<nw; ++iw)
        {
          PosType drOld =
            newpos[iw] - (W[iw]->R[iat] + oldScale[iw]*oldG[iw]);
          logGf_v[iw] = -m_oneover2tau * dot(drOld, drOld);
          rand_v[iw] = Random();
        }
        Psi.addRatio(W,iat,ratios,newG, newL);
        for(int iw=0; iw<nw; ++iw)
        {
          newScale[iw]   = getDriftScale(m_tauovermass,newG[iw]);
          PosType drNew  =
            (newpos[iw] + newScale[iw]*newG[iw]) - W[iw]->R[iat];
          RealType logGb =  -m_oneover2tau * dot(drNew, drNew);
          RealType x = logGb - logGf_v[iw];
          RealType prob = ratios[iw]*ratios[iw]*std::exp(x);
          if(rand_v[iw] < prob && ratios[iw] > 0.0)
          {
            accepted.push_back(W[iw]);
            nAccept++;
            W[iw]->R[iat] = newpos[iw];
            W[iw]->Age = 0;
            acc[iw] = true;
            R2acc[iw] += dot(delpos[iw], delpos[iw]);
          }
          else
            nReject++;
        }
        W.acceptMove_GPU(acc);
        if (accepted.size())
          Psi.update(accepted,iat);
      }
      for (int iw=0; iw < nw; iw++)
        if (W[iw]->Age)
          cerr << "Encountered stuck walker with iw=" << iw << endl;
      //	Psi.recompute(W, false);
      Psi.gradLapl(W, grad, lapl);
      H.evaluate (W, LocalEnergy, Txy);
      if (CurrentStep == 1)
        LocalEnergyOld = LocalEnergy;
      // Now, attempt nonlocal move
      accepted.clear();
      vector<int> iatList;
      vector<PosType> accPos;
      for (int iw=0; iw<nw; iw++)
      {
        int ibar = NLop.selectMove(Random(), Txy[iw]);
        // cerr << "Txy[iw].size() = " << Txy[iw].size() << endl;
        if (ibar)
        {
          accepted.push_back(W[iw]);
          int iat = Txy[iw][ibar].PID;
          iatList.push_back(iat);
          accPos.push_back(W[iw]->R[iat] + Txy[iw][ibar].Delta);
        }
      }
      if (accepted.size())
      {
        //   W.proposeMove_GPU(newpos, iatList);
        Psi.ratio(accepted,iatList, accPos, ratios, newG, newL);
        Psi.update(accepted,iatList);
        for (int i=0; i<accepted.size(); i++)
          accepted[i]->R[iatList[i]] = accPos[i];
        W.copyWalkersToGPU();
      }
      // Now branch
      for (int iw=0; iw<nw; iw++)
      {
        W[iw]->Weight *= branchEngine->branchWeight(LocalEnergy[iw], LocalEnergyOld[iw]);
        W[iw]->getPropertyBase()[R2ACCEPTED] = R2acc[iw];
        W[iw]->getPropertyBase()[R2PROPOSED] = R2prop[iw];
      }
      Mover->setMultiplicity(W.begin(), W.end());
      branchEngine->branch(CurrentStep,W);
      nw = W.getActiveWalkers();
      LocalEnergyOld.resize(nw);
      for (int iw=0; iw<nw; iw++)
        LocalEnergyOld[iw] = W[iw]->getPropertyBase()[LOCALENERGY];
    }
    while(step<nSteps);
    Psi.recompute(W, true);
    double accept_ratio = (double)nAccept/(double)(nAccept+nReject);
    Estimators->stopBlock(accept_ratio);
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    ++block;
    recordBlock(block);
  }
  while(block<nBlocks);
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
    Mover = new DMCUpdatePbyPWithRejection(W,Psi,H,Random);
    Mover->resetRun(branchEngine,Estimators);
    //Mover->initWalkersForPbyP(W.begin(),W.end());
  }
  else
  {
    int nw_multi=branchEngine->resetRun(qmcNode);
    if(nw_multi>1)
    {
      app_log() << " Current population " << W.getActiveWalkers() << " " << W.getGlobalNumWalkers()  << endl;
      app_log() << " The target population has changed. Multiply walkers by " << nw_multi << endl;
      W.createWalkers((nw_multi-1)*W.getActiveWalkers());
      setWalkerOffsets();
      app_log() << " New population " << W.getActiveWalkers() << " " << W.getGlobalNumWalkers()  << endl;
    }
  }

  branchEngine->checkParameters(W);

  //    Mover->updateWalkers(W.begin(),W.end());
  app_log() << "  DMC PbyP Update with a fluctuating population" << endl;
  Mover->MaxAge=1;
  app_log() << "  Steps per block = " << nSteps << endl;
  app_log() << "  Number of blocks = " << nBlocks << endl;
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
            << pool.getTotalSize() * sizeof(CudaRealType)
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
  vector<RealType> logPsi(W.WalkerList.size(), 0.0);
  //Psi.evaluateLog(W, logPsi);
  Psi.recompute(W, true);
  Estimators->start(nBlocks, true);
}

bool
DMCcuda::put(xmlNodePtr q)
{
  //nothing to add
  NLop.put(q);

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

/***************************************************************************
 * $RCSfile: DMCParticleByParticle.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: DMCcuda.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/

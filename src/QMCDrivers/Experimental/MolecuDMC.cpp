//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCDrivers/MolecuDMC.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

MolecuDMC::MolecuDMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h),
  KillNodeCrossing(0), KillWalker("no"),BranchInfo("default")
{
  RootName = "dmc";
  QMCType ="dmc";
  m_param.add(KillWalker,"killnode","string");
}

MolecuDMC::~MolecuDMC()
{
}

void MolecuDMC::setBranchInfo(const std::string& afile)
{
  BranchInfo=afile;
}

bool MolecuDMC::put(xmlNodePtr cur)
{
  return true;
}

/** Advance the walkers nblocks*nsteps timesteps.
 *
 * For each block:
 * <ul>
 *  <li> Advance walkers for nsteps
 *  For each timestep:
 *   <ul>
 *   <li> Move all the particles of a walker.
 *   <li> Calculate the properties for the new walker configuration.
 *   <li> Accept/reject the new configuration.
 *   <li> Accumulate the estimators.
 *   <li> Update the trial energy \f$ E_T \f$
 *   <li> Branch the population of walkers (birth/death algorithm).
 *   </ul>
 * <li> Flush the estimators and print to file.
 * <li> Update the estimate of the local energy.
 * <li> (Optional) Print the ensemble of walker configurations.
 * </ul>
 * Default mode: Print the ensemble of walker configurations
 * at the end of the run.
 */
bool MolecuDMC::run()
{
  KillNodeCrossing = (KillWalker == "yes");
  if(KillNodeCrossing)
  {
    app_log() << "Walkers will be killed if a node crossing is detected." << std::endl;
  }
  else
  {
    app_log() << "Walkers will be kept even if a node crossing is detected." << std::endl;
  }
  //add columns
  IndexType PopIndex = Estimators->addColumn("Population");
  IndexType EtrialIndex = Estimators->addColumn("Etrial");
  Estimators->reportHeader(AppendRun);
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  while(it != it_end)
  {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1;
    ++it;
  }
  m_oneover2tau = 0.5/Tau;
  m_sqrttau = sqrt(Tau);
  IndexType block = 0;
  int Population = W.getActiveWalkers();
  int tPopulation = W.getActiveWalkers();
  RealType Eest = branchEngine->E_T;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  do
  {
    IndexType step = 0;
    IndexType pop_acc=0;
    Estimators->startBlock();
    do
    {
      pop_acc += W.getActiveWalkers();
      if(KillNodeCrossing)
        advanceKillNodeCrossing(*branchEngine);
      else
        advanceRejectNodeCrossing(*branchEngine);
      step++;
      CurrentStep++;
      Estimators->accumulate(W);
      Eest = branchEngine->update(W.getActiveWalkers(), Eest);
      branchEngine->branch(CurrentStep,W);
    }
    while(step<nSteps);
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    Estimators->setColumn(PopIndex,static_cast<RealType>(pop_acc)/static_cast<RealType>(nSteps));
    Estimators->setColumn(EtrialIndex,Eest);
    Eest = Estimators->average(0);
    nAccept = 0;
    nReject = 0;
    block++;
    recordBlock(block);
    //create an output engine: could accumulate the configurations
    //if(block%Period4CheckPoint == 0) {
    //  HDFWalkerOutput WO(RootName,false,0);
    //  WO.get(W);
    //  WO.write(*branchEngine);
    //}
    W.reset();
  }
  while(block<nBlocks);
  //Need MPI-IO
  app_log() << "\t ratio = " << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot) << std::endl;
  Estimators->finalize();
  return true;
}

/**  Advance all the walkers one timstep while killing the walkers when node crossing is detected.
 * @param Branch class that controls the trial energy and branching
 */
template<class BRANCHER>
void
MolecuDMC::advanceKillNodeCrossing(BRANCHER& Branch)
{
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  while(it != it_end)
  {
    Walker_t& thisWalker(**it);
    thisWalker.Weight= 1.0;
    thisWalker.Multiplicity=1;
    //save old local energy
    RealType eold    = thisWalker.Properties(LOCALENERGY);
    RealType signold = thisWalker.Properties(SIGN);
    RealType emixed  = eold;
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandom(deltaR);
    W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
    //update the distance table associated with W
    //DistanceTable::update(W);
    W.update();
    //evaluate wave function
    ValueType logpsi(Psi.evaluateLog(W));
    bool accepted=false;
    if(Branch(Psi.getSign(),thisWalker.Properties(SIGN)))
    {
      thisWalker.Age++;
      thisWalker.willDie();
    }
    else
    {
      RealType enew(H.evaluate(W));
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      ValueType vsq = Dot(W.G,W.G);
      //converting gradients to drifts, D = tau*G (reuse G)
      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      drift = scale*W.G;
      deltaR = (*it)->R - W.R - drift;
      RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob= std::min(exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
      if(Random() > prob)
      {
        thisWalker.Age++;
      }
      else
      {
        accepted=true;
        thisWalker.R = W.R;
        thisWalker.Drift = drift;
        thisWalker.resetProperty(logpsi,Psi.getSign(),enew);
        H.saveProperty(thisWalker.getPropertyBase());
        emixed = (emixed+enew)*0.5;
        eold=enew;
      }
      //calculate the weight and multiplicity
      ValueType M = Branch.branchGF(Tau,emixed,0.0); //1.0-prob);
      if(thisWalker.Age > 3)
        M = std::min(0.5,M);
      else
        if(thisWalker.Age > 0)
          M = std::min(1.0,M);
      thisWalker.Weight = M;
      thisWalker.Multiplicity = M + Random();
    }
    Branch.accumulate(eold,thisWalker.Weight);
    if(accepted)
      ++nAccept;
    else
      ++nReject;
    ++it;
  }
}

/**  Advance all the walkers one timstep while killing the walkers when node crossing is detected.
 * @param Branch class that controls the trial energy and branching
 */
template<class BRANCHER>
void
MolecuDMC::advanceRejectNodeCrossing(BRANCHER& Branch)
{
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  while(it != it_end)
  {
    Walker_t& thisWalker(**it);
    thisWalker.Weight= 1.0;
    thisWalker.Multiplicity=1;
    //save old local energy
    RealType eold    = thisWalker.Properties(LOCALENERGY);
    RealType signold = thisWalker.Properties(SIGN);
    RealType emixed  = eold;
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandom(deltaR);
    W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
    //update the distance table associated with W
    //DistanceTable::update(W);
    W.update();
    //evaluate wave function
    ValueType logpsi(Psi.evaluateLog(W));
    bool accepted=false;
    if(Branch(Psi.getSign(),thisWalker.Properties(SIGN)))
    {
      thisWalker.Age++;
    }
    else
    {
      RealType enew(H.evaluate(W));
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      ValueType vsq = Dot(W.G,W.G);
      //converting gradients to drifts, D = tau*G (reuse G)
      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      drift = scale*W.G;
      deltaR = (*it)->R - W.R - drift;
      RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob= std::min(exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
      if(Random() > prob)
      {
        thisWalker.Age++;
      }
      else
      {
        accepted=true;
        thisWalker.R = W.R;
        thisWalker.Drift = drift;
        thisWalker.resetProperty(logpsi,Psi.getSign(),enew);
        H.saveProperty(thisWalker.getPropertyBase());
        emixed = (emixed+enew)*0.5;
        eold=enew;
      }
    }
    //calculate the weight and multiplicity
    ValueType M = Branch.branchGF(Tau,emixed,0.0);
    if(thisWalker.Age > 3)
      M = std::min(0.5,M);
    else
      if(thisWalker.Age > 0)
        M = std::min(1.0,M);
    thisWalker.Weight = M;
    thisWalker.Multiplicity = M + Random();
    Branch.accumulate(eold,thisWalker.Weight);
    if(accepted)
      ++nAccept;
    else
      ++nReject;
    ++it;
  }
}
}

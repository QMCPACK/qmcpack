//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMC/DMCPeta.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/Timer.h"

namespace qmcplusplus
{

/// Constructor.
DMCPeta::DMCPeta(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h), KillNodeCrossing(0), Reconfiguration("no"), Mover(0), BranchInterval(1)
{
  RootName = "dummy";
  QMCType ="DMCPeta";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_MULTIPLE,0);
  m_param.add(KillWalker,"killnode","string");
  m_param.add(BenchMarkRun,"benchmark","string");
  m_param.add(Reconfiguration,"reconfiguration","string");
  m_param.add(BranchInterval,"branchInterval","int");
  m_param.add(NonLocalMove,"nonlocalmove","string");
  m_param.add(NonLocalMove,"nonlocalmoves","string");
}

bool DMCPeta::run()
{
  resetUpdateEngine();
  //estimator is ready to collect data
  Estimators->setCollectionMode(true);
  Estimators->start(nBlocks,true);
  Timer myclock;
  IndexType block = 0;
  IndexType numPtcls=W.getTotalNum();
  RealType oneover2tau = 0.5/Tau;
  RealType sqrttau = std::sqrt(Tau);
  //temporary data to store differences
  ParticleSet::ParticlePos_t deltaR(numPtcls);
  ParticleSet::ParticleGradient_t G(numPtcls), dG(numPtcls);
  ParticleSet::ParticleLaplacian_t L(numPtcls), dL(numPtcls);
  CurrentStep = 0;
  typedef MCWalkerConfiguration::iterator WalkerIter_t;
  do // block
  {
    Mover->startBlock(nSteps);
    for(IndexType step=0; step< nSteps; step++, CurrentStep+=BranchInterval)
    {
      for(IndexType interval=0; interval<BranchInterval; ++interval)
      {
        for(WalkerIter_t it=W.begin(); it != W.end(); ++it)
        {
          //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
          Walker_t& thisWalker(**it);
          Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
          W.R = thisWalker.R;
          w_buffer.rewind();
          W.copyFromBuffer(w_buffer);
          Psi.copyFromBuffer(W,w_buffer);
          //create a 3N-Dimensional Gaussian with variance=1
          makeGaussRandomWithEngine(deltaR,Random);
          int nAcceptTemp(0);
          int nRejectTemp(0);
          RealType eold(thisWalker.Properties(LOCALENERGY));
          RealType enew(eold);
          RealType rr_proposed=0.0;
          RealType rr_accepted=0.0;
          //loop over particles
          for(int iat=0; iat<numPtcls; iat++)
          {
            PosType dr(sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
            PosType newpos(W.makeMove(iat,dr));
            RealType ratio=Psi.ratio(W,iat,dG,dL);
            RealType rr=dot(dr,dr);
            rr_proposed+=rr;
            if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
            {
              //node crossing detected
              ++nRejectTemp;
              W.rejectMove(iat);
              Psi.rejectMove(iat);
            }
            else
            {
              G = W.G+dG;
              RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
              RealType scale=getDriftScale(Tau,G);
              dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);
              RealType logGb = -oneover2tau*dot(dr,dr);
              RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
              if(Random() < prob)
              {
                //move is accepted
                ++nAcceptTemp;
                W.acceptMove(iat);
                Psi.acceptMove(W,iat);
                W.G = G;
                W.L += dL;
                assignDrift(scale,G,thisWalker.Drift);
                rr_accepted+=rr;
              }
              else
              {
                ++nRejectTemp;
                W.rejectMove(iat);
                Psi.rejectMove(iat);
              }
            }
          }//for(int iat=0; iat<NumPtcl; iat++)
          if(nAcceptTemp>0)
          {
            //need to overwrite the walker properties
            thisWalker.R = W.R;
            w_buffer.rewind();
            W.copyToBuffer(w_buffer);
            //RealType psi = Psi.evaluate(W,w_buffer);
            RealType logpsi = Psi.evaluateLog(W,w_buffer);
            enew= H.evaluate(W);
            //thisWalker.resetProperty(std::log(abs(psi)),psi,enew,rr_accepted,rr_proposed,1.0);
            thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,1.0 );
            H.auxHevaluate(W,thisWalker);
            H.saveProperty(thisWalker.getPropertyBase());
          }
          else
          {
            thisWalker.rejectedMove();
            thisWalker.Age++;
            rr_accepted=0.0;
            enew=eold;//copy back old energy
          }
          thisWalker.Weight *= branchEngine->branchWeight(eold,enew);
          nAccept += nAcceptTemp;
          nReject += nRejectTemp;
        }//for(WalkerIter_t it=W.begin();it != W.end(); ++it)
      }//interval
      //calculate multiplicity based on the weights: see QMCUpdateBase::setMultiplicity
      Mover->setMultiplicity(W.begin(),W.end());
      //time to branch: see SimpleFixedNodeBranch::branch
      branchEngine->branch(CurrentStep,W);
      if(storeConfigs && (CurrentStep%storeConfigs == 0))
      {
        ForwardWalkingHistory.storeConfigsForForwardWalking(W);
        W.resetWalkerParents();
      }
    }//steps
    ++block;
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    recordBlock(block);
  }
  while(block<nBlocks &&  myclock.elapsed()<MaxCPUSecs);
  Estimators->stop();
  return finalize(block);
}

void DMCPeta::resetUpdateEngine()
{
  bool fixW=false;
  if(Mover==0) //disable switching update modes for DMCPeta in a run
  {
    branchEngine->initWalkerController(W,Tau,fixW);
    Mover= new DMCUpdatePbyPWithRejection(W,Psi,H,Random);
    Mover->MaxAge=3;
    Mover->resetRun(branchEngine,Estimators);
    Mover->initWalkersForPbyP(W.begin(),W.end());
  }
  app_log() << "  BranchInterval = " << BranchInterval << endl;
  app_log() << "  Steps per block = " << nSteps << endl;
  app_log() << "  Number of blocks = " << nBlocks << endl;
}

bool DMCPeta::put(xmlNodePtr q)
{
  return true;
}

}

/***************************************************************************
 * $RCSfile: DMC.cpp,v $   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: DMC.cpp 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/

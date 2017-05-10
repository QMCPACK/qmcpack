//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
//#define TEST_INNERBRANCH


namespace qmcplusplus
{

/// Constructor.
RNDMCUpdatePbyPCeperley::RNDMCUpdatePbyPCeperley(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,h,rg), maxS(100), estimateCrossings(0)
{
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::advance")); //timer for the walker loop
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::movePbyP")); //timer for MC, ratio etc
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::updateMBO")); //timer for measurements
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::energy")); //timer for measurements
  TimerManager.addTimer(myTimers[0]);
  TimerManager.addTimer(myTimers[1]);
  TimerManager.addTimer(myTimers[2]);
  TimerManager.addTimer(myTimers[3]);
  myParams.add(maxS,"Smax","int");
  myParams.add(estimateCrossings,"estcross","int");
  myParams.add(efn,"fnenergy","float");
}

/// destructor
RNDMCUpdatePbyPCeperley::~RNDMCUpdatePbyPCeperley() { }


void RNDMCUpdatePbyPCeperley::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=true;
  for (; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,UpdatePbyP);
    Walker_t::Buffer_t tbuffer;
    RealType logpsi=Psi.registerData(W,tbuffer);
    thisWalker.DataSet=tbuffer;
    //setScaledDriftPbyP(m_tauovermass,W.G,(*it)->Drift);
//         RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    RealType bene = H.evaluate(W);
    thisWalker.resetProperty(logpsi,Psi.getPhase(),bene, 0.0,0.0, 1.0);
    H.saveProperty(thisWalker.getPropertyBase());
    ValueType altR = Psi.alternateRatio(W);
    if (altR<0)
    {
      PosType p1(W.R[0]);
      W.R[0]=W.R[1];
      W.R[1]=p1;
      logpsi=Psi.updateBuffer(W,thisWalker.DataSet,true);
      W.saveWalker(thisWalker);
      altR = Psi.alternateRatio(W);
    }
    RealType fene  = -0.5*(Sum(W.L)+Dot(W.G,W.G)) + thisWalker.Properties(LOCALPOTENTIAL);
    thisWalker.resetReleasedNodeProperty(bene,fene,altR);
    thisWalker.ReleasedNodeAge=0;
    thisWalker.ReleasedNodeWeight=0.0;
  }
}


/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void RNDMCUpdatePbyPCeperley::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end
    , bool measure)
{
  myTimers[0]->start();
  for(; it != it_end; ++it)
  {
    //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
    Walker_t& thisWalker(**it);
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    W.loadWalker(thisWalker,true);
    //W.R = thisWalker.R;
    //w_buffer.rewind();
    //W.copyFromBuffer(w_buffer);
    Psi.copyFromBuffer(W,w_buffer);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    int nAcceptTemp(0);
    int nRejectTemp(0);
    //copy the old energy and scale factor of drift
    RealType feold(thisWalker.Properties(ALTERNATEENERGY));
    RealType beold(thisWalker.Properties(LOCALENERGY));
    RealType oldAltR(thisWalker.Properties(SIGN));
    RealType bene(beold);
    RealType rr_proposed=0.0;
    RealType rr_accepted=0.0;
    RealType gf_acc=1.0;
    nNodeCrossing=0;
    bool inFNpop(thisWalker.ReleasedNodeAge==0);
    myTimers[1]->start();
    for(int iat=0; iat<NumPtcl; ++iat)
    {
      bool crossattempt(false);
      //get the displacement
      GradType grad_iat=Psi.evalGrad(W,iat);
      ParticleSet::ParticleGradient_t pg1(W.G);
      pg1=0;
      if (inFNpop)
        Psi.alternateGrad(pg1);
      PosType dr;
      dr = m_tauovermass*grad_iat+m_sqrttau*deltaR[iat];
      //RealType rr=dot(dr,dr);
      RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
      rr_proposed+=rr;
      if(rr>m_r2max)
      {
        ++nRejectTemp;
        continue;
      }
      //PosType newpos(W.makeMove(iat,dr));
      if(!W.makeMoveAndCheck(iat,dr))
        continue;
      PosType newpos(W.R[iat]);
      RealType ratio = Psi.ratioGrad(W,iat,grad_iat);
      bool valid_move=false;
      if (branchEngine->phaseChanged(Psi.getAlternatePhaseDiff(iat)))
      {
        crossattempt=true;
      }
      else
        if ((inFNpop)&&(estimateCrossings))
        {
          ParticleSet::ParticleGradient_t pg2(W.G);
          pg2=0;
          Psi.alternateGrad(pg2);
          if (dot(pg2[iat],pg1[iat])>0)
          {
            RealType graddot=dot(pg1[iat],pg1[iat]);
            RealType graddot2=dot(pg2[iat],pg2[iat]);
            RealType magG=std::sqrt(graddot*graddot2);
            graddot = -4.0*m_oneover2tau/magG;
            RealType crossing=1.0-std::exp(graddot);
            if (RandomGen()>crossing)
              crossattempt=true;
          }
        }
      RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
//           getScaledDrift(m_tauovermass, grad_iat, dr);
      dr = thisWalker.R[iat] - newpos - m_tauovermass*grad_iat;
      RealType logGb = -m_oneover2tau*dot(dr,dr);
      RealType prob = ratio*ratio*std::exp(logGb-logGf);
      if(RandomGen() < prob)
      {
        valid_move=true;
        ++nAcceptTemp;
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
        rr_accepted+=rr;
        gf_acc *=prob;//accumulate the ratio
        if (crossattempt)
          ++nNodeCrossing;
      }
      else
      {
        ++nRejectTemp;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
      }
    }
    myTimers[1]->stop();
    //RealType nodecorr_old=thisWalker.Properties(DRIFTSCALE);
    //RealType nodecorr=nodecorr_old;
    bool advanced=true;
    if(nAcceptTemp>0)
    {
      //need to overwrite the walker properties
      myTimers[2]->start();
      thisWalker.Age++;
      thisWalker.R = W.R;
      //w_buffer.rewind();
      //W.updateBuffer(w_buffer);
      //         RealType logpsi = Psi.evaluateLog(W,w_buffer);
//         RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
      RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
      W.saveWalker(thisWalker);
      myTimers[2]->stop();
      myTimers[3]->start();
      bene = H.evaluate(W);
      myTimers[3]->stop();
      //nodecorr=getNodeCorrection(W.G,thisWalker.Drift);
      //thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),bene,rr_accepted,rr_proposed,1.0 );
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    else
    {
      //all moves are rejected: does not happen normally with reasonable wavefunctions
      advanced=false;
      thisWalker.Age++;
      thisWalker.Properties(R2ACCEPTED)=0.0;
      H.rejectedMove(W,thisWalker);
      ++nAllRejected;
      bene =beold;//copy back old energy
      gf_acc=1.0;
    }
//       RealType sigma_last = thisWalker.ReleasedNodeWeight/std::abs(thisWalker.ReleasedNodeWeight);
    ValueType altR = Psi.alternateRatio(W);
    RealType fenew = -0.5*(Sum(W.L)+Dot(W.G,W.G)) + thisWalker.Properties(LOCALPOTENTIAL);
    thisWalker.resetReleasedNodeProperty(bene ,fenew,altR);
    RealType bareBranchWeight=branchEngine->branchWeightBare(bene,beold);
    RealType fixedBranchWeight=branchEngine->branchWeightReleasedNode(bene,beold,efn);
    if ((nNodeCrossing==0)&&(thisWalker.ReleasedNodeAge==0))
    {
      thisWalker.ReleasedNodeWeight=altR;
//         thisWalker.Weight *= branchEngine->branchWeightReleasedNode(bene,beold,efn);
      thisWalker.Weight *= bareBranchWeight;
    }
    else
      if (thisWalker.ReleasedNodeAge<maxS)
      {
        //Case where walker is in released node population
//         thisWalker.Weight *= branchEngine->branchWeightReleasedNode(enew,eold,efn);
        thisWalker.Weight *= fixedBranchWeight;
//         thisWalker.Weight *= bareBranchWeight;
        thisWalker.ReleasedNodeAge++;
        thisWalker.ReleasedNodeWeight=altR;
      }
      else
      {
//         case where walker has exceeded the time limit.
        thisWalker.Weight=0;
        thisWalker.ReleasedNodeWeight=0;
        thisWalker.ReleasedNodeAge++;
      }
    nAccept += nAcceptTemp;
    nReject += nRejectTemp;
  }
  myTimers[0]->stop();
}

}


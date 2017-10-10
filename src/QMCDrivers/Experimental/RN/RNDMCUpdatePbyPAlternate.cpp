//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
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
RNDMCUpdatePbyPAlternate::RNDMCUpdatePbyPAlternate(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,h,rg), maxS(100), estimateCrossings(0), maxcopy(10)
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
  myParams.add(maxcopy,"maxCopy","int");
  myParams.add(efn,"fnenergy","float");
  myParams.add(estimateCrossings,"estcross","int");
}

/// destructor
RNDMCUpdatePbyPAlternate::~RNDMCUpdatePbyPAlternate() { }

void RNDMCUpdatePbyPAlternate::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
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
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    RealType fene = H.evaluate(W);
    thisWalker.resetProperty(logpsi,0.0,fene, 0.0,0.0, nodecorr);
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
    thisWalker.ReleasedNodeAge=0;
//         thisWalker.Weight=altR;
    RealType bene = -0.5*(Sum(W.L)+Dot(W.G,W.G)) + thisWalker.Properties(LOCALPOTENTIAL);
    thisWalker.resetReleasedNodeProperty(fene ,bene ,altR);
  }
}

/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void RNDMCUpdatePbyPAlternate::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end
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
    RealType beold(thisWalker.Properties(ALTERNATEENERGY));
    RealType eold(thisWalker.Properties(LOCALENERGY));
    RealType oldAltR(thisWalker.Properties(SIGN));
    RealType vqold(thisWalker.Properties(DRIFTSCALE));
    RealType enew(eold);
    RealType rr_proposed=0.0;
    RealType rr_accepted=0.0;
    RealType gf_acc=1.0;
    nNodeCrossing=0;
    bool inFNpop = (thisWalker.ReleasedNodeAge==0);
    myTimers[1]->start();
    for(int iat=0; iat<NumPtcl; ++iat)
    {
      bool attemptnodecrossing(false);
      //get the displacement
      GradType grad_iat;
//         if (inFNpop)
      grad_iat=Psi.evalGrad(W,iat);
//         else
//           grad_iat=Psi.alternateEvalGrad(W,iat);
      PosType dr;
      //RealType sc=getDriftScale(m_tauovermass,grad_iat);
      //PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_iat));
//         getScaledDrift(m_tauovermass, grad_iat, dr);
//         dr += m_sqrttau * deltaR[iat];
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
      RealType ratio;
//         GradType grad_iat2(grad_iat);
//         if (inFNpop)
      ratio = Psi.ratioGrad(W,iat,grad_iat);
//         else
//           ratio = Psi.alternateRatioGrad(W,iat,grad_iat);
      bool valid_move=false;
      //node is crossed reject the move
      //if(Psi.getPhase() > std::numeric_limits<RealType>::epsilon())
      //if(branchEngine->phaseChanged(Psi.getPhase(),thisWalker.Properties(SIGN)))
//         if (inFNpop)
      if (branchEngine->phaseChanged(Psi.getAlternatePhaseDiff(iat)))
        attemptnodecrossing=true;
//         else if ((inFNpop)&&(estimateCrossings))
//         {
//           if (dot(grad_iat2,grad_iat)>0)
//           {
//             RealType graddot=dot(grad_iat,grad_iat);
//             RealType graddot2=dot(grad_iat2,grad_iat2);
//             RealType magG=std::sqrt(graddot*graddot2);
// //           if (graddot>0)
// //           {
//             graddot = -4.0*m_oneover2tau/magG;
//             RealType crossing=1.0-std::exp(graddot);
//             if (RandomGen()>crossing) attemptnodecrossing=true;
//           }
//         }
      RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
      //Use the force of the particle iat
      //RealType scale=getDriftScale(m_tauovermass,grad_iat);
      //dr = thisWalker.R[iat]-newpos-scale*real(grad_iat);
//           getScaledDrift(m_tauovermass, grad_iat, dr);
//           dr = thisWalker.R[iat] - newpos - dr;
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
        if (attemptnodecrossing)
          nNodeCrossing++;
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
      RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
      W.saveWalker(thisWalker);
      myTimers[2]->stop();
      myTimers[3]->start();
      enew= H.evaluate(W);
      myTimers[3]->stop();
      //nodecorr=getNodeCorrection(W.G,thisWalker.Drift);
      //thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,1.0 );
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
      enew=eold;//copy back old energy
      gf_acc=1.0;
    }
//       RealType sigma_last = thisWalker.ReleasedNodeWeight/std::abs(thisWalker.ReleasedNodeWeight);
    ValueType altR = Psi.alternateRatio(W);
    RealType benew = -0.5*(Sum(W.L)+Dot(W.G,W.G)) + thisWalker.Properties(LOCALPOTENTIAL);
    thisWalker.resetReleasedNodeProperty(enew,benew,altR);
    RealType rnweight=branchEngine->branchWeightReleasedNode(benew,beold,efn);
    RealType fnweight=branchEngine->branchWeight(enew,eold);
    if ((thisWalker.ReleasedNodeAge==0)&&(nNodeCrossing==0))
    {
//         Case where walker is still in fixed node population
      thisWalker.Weight *= (altR/oldAltR)*fnweight;
      thisWalker.ReleasedNodeWeight = 1.0/altR;
    }
    else
      if ( thisWalker.ReleasedNodeAge==0)
      {
        //First crossing.
        thisWalker.Weight *= altR/oldAltR*(fnweight);
        thisWalker.ReleasedNodeWeight = 1.0/altR;
        thisWalker.ReleasedNodeAge++;
      }
      else
        if ( thisWalker.ReleasedNodeAge<maxS )
        {
          //Case where walker is in released node population
          thisWalker.Weight *= rnweight;
          thisWalker.ReleasedNodeWeight = 1.0/altR;
          thisWalker.ReleasedNodeAge++;
        }
        else
        {
          //         case where walker has exceeded the time limit.
          thisWalker.Weight = 0;
          thisWalker.ReleasedNodeWeight = 0;
          thisWalker.ReleasedNodeAge++;
        }
    nAccept += nAcceptTemp;
    nReject += nRejectTemp;
  }
  myTimers[0]->stop();
}

}


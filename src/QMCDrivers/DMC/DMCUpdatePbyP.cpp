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
DMCUpdatePbyPWithRejection::DMCUpdatePbyPWithRejection(MCWalkerConfiguration& w
    , TrialWaveFunction& psi, QMCHamiltonian& h , RandomGenerator_t& rg)
  : QMCUpdateBase(w,psi,h,rg)
{
  myTimers.push_back(new NewTimer("DMCUpdatePbyP::advance")); //timer for the walker loop
  myTimers.push_back(new NewTimer("DMCUpdatePbyP::movePbyP")); //timer for MC, ratio etc
  myTimers.push_back(new NewTimer("DMCUpdatePbyP::updateMBO")); //timer for measurements
  myTimers.push_back(new NewTimer("DMCUpdatePbyP::energy")); //timer for measurements
  TimerManager.addTimer(myTimers[0]);
  TimerManager.addTimer(myTimers[1]);
  TimerManager.addTimer(myTimers[2]);
  TimerManager.addTimer(myTimers[3]);
}

/// destructor
DMCUpdatePbyPWithRejection::~DMCUpdatePbyPWithRejection() { }

/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void DMCUpdatePbyPWithRejection::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
    bool measure)
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
    RealType eold(thisWalker.Properties(LOCALENERGY));
    RealType vqold(thisWalker.Properties(DRIFTSCALE));
    RealType enew(eold);
    RealType rr_proposed=0.0;
    RealType rr_accepted=0.0;
    RealType gf_acc=1.0;
    myTimers[1]->start();
    for(int iat=0; iat<NumPtcl; ++iat)
    {
      PosType dr;
      //get the displacement
      //RealType sc=getDriftScale(m_tauovermass,W.G[iat]);
      //PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));
      getScaledDrift(m_tauovermass,W.G[iat],dr);
      dr += m_sqrttau*deltaR[iat];
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
      RealType ratio=Psi.ratio(W,iat,dG,dL);
      bool valid_move=false;
      //node is crossed reject the move
      //if(Psi.getPhase() > numeric_limits<RealType>::epsilon())
      if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
      {
        ++nRejectTemp;
        ++nNodeCrossing;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
      }
      else
      {
        G = W.G+dG;
        RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
        //Use the force of the particle iat
        //RealType scale=getDriftScale(m_tauovermass,G[iat]);
        //dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);
        getScaledDrift(m_tauovermass, G[iat], dr);
        dr = thisWalker.R[iat] - newpos - dr;
        RealType logGb = -m_oneover2tau*dot(dr,dr);
        RealType prob = ratio*ratio*std::exp(logGb-logGf);
        //this is useless
        //RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
        if(RandomGen() < prob)
        {
          valid_move=true;
          ++nAcceptTemp;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
          W.G = G;
          W.L += dL;
          rr_accepted+=rr;
          gf_acc *=prob;//accumulate the ratio
        }
        else
        {
          ++nRejectTemp;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
    }
    myTimers[1]->stop();
    RealType nodecorr_old=thisWalker.Properties(DRIFTSCALE);
    RealType nodecorr=nodecorr_old;
    bool advanced=true;
    if(UseTMove)
      nonLocalOps.reset();
    if(nAcceptTemp>0)
    {
      //need to overwrite the walker properties
      myTimers[2]->start();
      thisWalker.Age=0;
      thisWalker.R = W.R;
      nodecorr=getNodeCorrection(W.G,drift);
      //w_buffer.rewind();
      //W.copyToBuffer(w_buffer);
      RealType logpsi = Psi.evaluateLog(W,w_buffer);
      W.saveWalker(thisWalker);
      myTimers[2]->stop();
      myTimers[3]->start();
      if(UseTMove)
        enew= H.evaluate(W,nonLocalOps.Txy);
      else
        enew= H.evaluate(W);
      myTimers[3]->stop();
      //thisWalker.resetProperty(std::log(abs(psi)),psi,enew,rr_accepted,rr_proposed,nodecorr);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr );
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    else
    {
      //all moves are rejected: does not happen normally with reasonable wavefunctions
      advanced=false;
      H.rejectedMove(W,thisWalker);
      thisWalker.Age++;
      thisWalker.Properties(R2ACCEPTED)=0.0;
      ++nAllRejected;
      enew=eold;//copy back old energy
      gf_acc=1.0;
    }
    if(UseTMove)
    {
      //make a non-local move
      int ibar=nonLocalOps.selectMove(RandomGen());
      if(ibar)
      {
        myTimers[2]->start();
        int iat=nonLocalOps.id(ibar);
        PosType newpos(W.makeMove(iat, nonLocalOps.delta(ibar)));
        RealType ratio=Psi.ratio(W,iat,dG,dL);
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
        W.G += dG;
        W.L += dL;
        //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
        //thisWalker.R[iat]=W.R[iat];
        //w_buffer.rewind();
        //W.copyToBuffer(w_buffer);
        RealType logpsi = Psi.evaluateLog(W,w_buffer);
        W.saveWalker(thisWalker);
        ++NonLocalMoveAccepted;
        myTimers[2]->stop();
      }
    }
    //2008-06-26: select any
    //bare green function by setting nodecorr=nodecorr_old=1.0
    thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
    //Filtering extreme energies
    //thisWalker.Weight *= branchEngine->branchWeight(eold,enew);
    //using the corrections: see QMCUpdateBase::getNodeCorrection
    //thisWalker.Weight *= branchEngine->branchWeight(enew,eold,nodecorr,nodecorr_old);
    //using the corrections: see QMCUpdateBase::getNodeCorrection  including gf_acc
    //RealType odd=std::min(gf_acc,1.0)
    //thisWalker.Weight *= branchEngine->branchWeight(enew,eold,nodecorr,nodecorr_oldi,odd);
    nAccept += nAcceptTemp;
    nReject += nRejectTemp;
  }
  myTimers[0]->stop();
}

/// Constructor.
DMCUpdatePbyPWithKill::DMCUpdatePbyPWithKill(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
    RandomGenerator_t& rg): QMCUpdateBase(w,psi,h,rg)
{ }

/// destructor
DMCUpdatePbyPWithKill::~DMCUpdatePbyPWithKill() { }

/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void DMCUpdatePbyPWithKill::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
    bool measure)
{
  app_error() << "  DMCUpdatePbyPWithKill::advanceWalkers in not implemented." << endl;
}

}

/***************************************************************************
 * $RCSfile: DMCUpdatePbyP.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

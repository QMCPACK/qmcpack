//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
#include "QMCDrivers/DriftOperators.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus { 

  /// Constructor.
  DMCNonLocalUpdate::DMCNonLocalUpdate(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      RandomGenerator_t& rg): QMCUpdateBase(w,psi,h,rg)
    { }
  
  /// destructor
  DMCNonLocalUpdate::~DMCNonLocalUpdate() { }

  /** advance all the walkers with killnode==no
   * @param nat number of particles to move
   * 
   * When killnode==no, any move resulting in node-crossing is treated
   * as a normal rejection.
   */
  void DMCNonLocalUpdate::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
      bool measure) 
  {

    //RealType plusFactor(Tau*Gamma);
    //RealType minusFactor(-Tau*(1.0-Alpha*(1.0+Gamma)));

    while(it != it_end) {
      Walker_t& thisWalker(**it);

      //save old local energy
      RealType eold    = thisWalker.Properties(LOCALENERGY);
      RealType signold = thisWalker.Properties(SIGN);
      RealType enew  = eold;

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(deltaR,RandomGen);

      W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
      
      //update the distance table associated with W
      //DistanceTable::update(W);
      W.update();
      
      //evaluate wave function
      RealType logpsi(Psi.evaluateLog(W));

      nonLocalOps.reset();

      bool accepted=false; 
      if(branchEngine->phaseChanged(Psi.getPhase(),thisWalker.Properties(SIGN))) {
        thisWalker.Age++;
        ++nReject;
      } else {
        //RealType enew(H.evaluate(W,nonLocalOps.Txy));
        enew=H.evaluate(W,nonLocalOps.Txy);
        RealType logGf = -0.5*Dot(deltaR,deltaR);
        setScaledDrift(Tau,W.G,drift);

        deltaR = (*it)->R - W.R - drift;
        RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);

        RealType prob= std::min(std::exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
        if(RandomGen() > prob){
          thisWalker.Age++;
          ++nReject;
          enew=eold;
        } else {
          accepted=true;  
          thisWalker.R = W.R;
          thisWalker.Drift = drift;
          thisWalker.resetProperty(logpsi,Psi.getPhase(),enew);
          H.saveProperty(thisWalker.getPropertyBase());
          //emixed = (emixed+enew)*0.5;
          //eold=enew;
          ++nAccept;
        }
      }

      int ibar=nonLocalOps.selectMove(RandomGen());

      //make a non-local move
      if(ibar) {
        int iat=nonLocalOps.id(ibar);
        W.R[iat] += nonLocalOps.delta(ibar);
        W.update();
        logpsi=Psi.evaluateLog(W);
        setScaledDrift(Tau,W.G,thisWalker.Drift);
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eold);
        thisWalker.R[iat] = W.R[iat];
        ++NonLocalMoveAccepted;
      } 

      //thisWalker.Weight *= branchEngine->branchGF(Tau,emixed,0.0);
      thisWalker.Weight *= branchEngine->branchWeight(Tau,eold,enew);
      //branchEngine->accumulate(eold,1);
      ++it;
    }
  }

  /// Constructor.
  DMCNonLocalUpdatePbyP::DMCNonLocalUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      RandomGenerator_t& rg): QMCUpdateBase(w,psi,h,rg)
    { }
  
  /// destructor
  DMCNonLocalUpdatePbyP::~DMCNonLocalUpdatePbyP() { }

  /** advance all the walkers with killnode==no
   * @param nat number of particles to move
   * 
   * When killnode==no, any move resulting in node-crossing is treated
   * as a normal rejection.
   */
  void DMCNonLocalUpdatePbyP::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
      bool measure) 
  {

    while(it != it_end) {
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      RealType eold(thisWalker.Properties(LOCALENERGY));
      //RealType emixed(eold), enew(eold);
      RealType enew(eold);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      bool notcrossed(true);
      int nAcceptTemp(0);
      int nRejectTemp(0);
      int iat=0;

      RealType rr_proposed=0.0;
      RealType rr_accepted=0.0;
      while(iat<NumPtcl) {//particle-by-particle move
        PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        PosType newpos(W.makeMove(iat,dr));

        RealType ratio=Psi.ratio(W,iat,dG,dL);

        RealType rr=dot(dr,dr);
        rr_proposed+=rr;

        //if(ratio < 0.0) {//node is crossed reject the move
        if(Psi.getPhase() > numeric_limits<RealType>::epsilon()) {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.rejectMove(iat); Psi.rejectMove(iat);
        } else {
          G = W.G+dG;
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
          RealType scale=getDriftScale(Tau,G);
          dr = thisWalker.R[iat]-newpos-scale*real(G[iat]); 
          RealType logGb = -m_oneover2tau*dot(dr,dr);
          RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
          if(RandomGen() < prob) { 
            ++nAcceptTemp;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
            W.G = G;
            W.L += dL;
            assignDrift(scale,G,thisWalker.Drift);
            rr_accepted+=rr;
          } else {
            ++nRejectTemp; 
            W.rejectMove(iat); Psi.rejectMove(iat);
          }
        } 
        ++iat;
      }//end of drift+diffusion

      nonLocalOps.reset();
      if(nAcceptTemp>0) {//need to overwrite the walker properties
        thisWalker.R = W.R;
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        RealType psi = Psi.evaluate(W,w_buffer);
        enew= H.evaluate(W,nonLocalOps.Txy);
        thisWalker.resetProperty(std::log(abs(psi)),psi,enew);
        H.saveProperty(thisWalker.getPropertyBase());
        //emixed = (eold+enew)*0.5e0;
      } else {
        thisWalker.Age++;
        ++nAllRejected;
        rr_accepted=0.0;
        enew=eold;//copy back old energy
      }

      int ibar = nonLocalOps.selectMove(RandomGen());

      //make a non-local move
      if(ibar) {
        int iat=nonLocalOps.id(ibar);
        PosType newpos(W.makeMove(iat, nonLocalOps.delta(ibar)));
        RealType ratio=Psi.ratio(W,iat,dG,dL);
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
        W.G += dG;
        W.L += dL;
        setScaledDrift(Tau,W.G,thisWalker.Drift);

        thisWalker.R[iat]=W.R[iat];
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        RealType psi = Psi.evaluate(W,w_buffer);

        ++NonLocalMoveAccepted;
      } 

      //thisWalker.Weight *= branchEngine->branchGF(Tau*rr_accepted/rr_proposed,emixed,0.0);
      thisWalker.Weight *= branchEngine->branchWeight(Tau*rr_accepted/rr_proposed,eold,enew);
      //branchEngine->accumulate(eold,1);
      nAccept += nAcceptTemp;
      nReject += nRejectTemp;
      ++it;
    }
  }
}

/***************************************************************************
 * $RCSfile: DMCNonLocalUpdate.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

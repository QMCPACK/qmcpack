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

    for(; it!=it_end; ++it)
    {
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
      if(branchEngine->phaseChanged(Psi.getPhase(),thisWalker.Properties(SIGN))) 
      {
        thisWalker.rejectedMove();
        thisWalker.Age++;
        ++nReject;
      } 
      else 
      {
        //RealType enew(H.evaluate(W,nonLocalOps.Txy));
        enew=H.evaluate(W,nonLocalOps.Txy);
        RealType logGf = -0.5*Dot(deltaR,deltaR);
        setScaledDrift(Tau,W.G,drift);

        deltaR = (*it)->R - W.R - drift;
        RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);

        RealType prob= std::min(std::exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
        if(RandomGen() > prob){
          thisWalker.Age++;
          thisWalker.rejectedMove();
          ++nReject;
          enew=eold;
        } else {
          accepted=true;  
          thisWalker.R = W.R;
          thisWalker.Drift = drift;
          thisWalker.resetProperty(logpsi,Psi.getPhase(),enew);
          H.auxHevaluate(W,thisWalker);
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

      thisWalker.Weight *= branchEngine->branchWeight(eold,enew);
      //branchEngine->accumulate(eold,1);
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

    for(; it!=it_end; ++it)
    {
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType vqold(thisWalker.Properties(DRIFTSCALE));
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

      RealType rr_proposed=0.0;
      RealType rr_accepted=0.0;
      for(int iat=0; iat<NumPtcl; ++iat)
      {
        //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        RealType sc=getDriftScale(Tau,W.G[iat]);
        PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));

        //RealType rr=dot(dr,dr);
        RealType rr=Tau*dot(deltaR[iat],deltaR[iat]);
        rr_proposed+=rr;

        if(rr>m_r2max)//too big move
        {
          ++nRejectTemp; continue;
        }

        PosType newpos(W.makeMove(iat,dr));
        RealType ratio=Psi.ratio(W,iat,dG,dL);

        //node is crossed reject the move
        if(Psi.getPhase() > numeric_limits<RealType>::epsilon()) 
        {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.rejectMove(iat); Psi.rejectMove(iat);
        } 
        else 
        {
          G = W.G+dG;
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);

          //RealType scale=getDriftScale(Tau,G);
          RealType scale=getDriftScale(Tau,G[iat]);
          dr = thisWalker.R[iat]-newpos-scale*real(G[iat]); 

          RealType logGb = -m_oneover2tau*dot(dr,dr);
          RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
          if(RandomGen() < prob) 
          { 
            ++nAcceptTemp;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
            W.G = G;
            W.L += dL;
            //assignDrift(scale,G,thisWalker.Drift);
            rr_accepted+=rr;
          } 
          else 
          {
            ++nRejectTemp; 
            W.rejectMove(iat); Psi.rejectMove(iat);
          }
        } 
      }//end of drift+diffusion for all the particles of a walker

      nonLocalOps.reset();
      if(nAcceptTemp>0) 
      {//need to overwrite the walker properties
        thisWalker.R = W.R;
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        //RealType psi = Psi.evaluate(W,w_buffer);
        RealType logpsi = Psi.evaluateLog(W,w_buffer);
        enew= H.evaluate(W,nonLocalOps.Txy);
        //thisWalker.resetProperty(std::log(abs(psi)),psi,enew,rr_accepted,rr_proposed,1.0);
        thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,1.0);
        H.auxHevaluate(W,thisWalker);
        H.saveProperty(thisWalker.getPropertyBase());

        thisWalker.Drift=W.G;
      } 
      else 
      {
        thisWalker.rejectedMove();
        thisWalker.Age++;
        thisWalker.Properties(R2ACCEPTED)=0.0;
        ++nAllRejected;
        enew=eold;//copy back old energy
      }

      int ibar = nonLocalOps.selectMove(RandomGen());

      //make a non-local move
      if(ibar) 
      {
        int iat=nonLocalOps.id(ibar);
        PosType newpos(W.makeMove(iat, nonLocalOps.delta(ibar)));
        RealType ratio=Psi.ratio(W,iat,dG,dL);
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
        W.G += dG;
        W.L += dL;

        PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);

        thisWalker.R[iat]=W.R[iat];
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        //RealType psi = Psi.evaluate(W,w_buffer);
        RealType logpsi = Psi.evaluateLog(W,w_buffer);

        ++NonLocalMoveAccepted;
      }

      thisWalker.Weight *= branchEngine->branchWeight(eold,enew);

      nAccept += nAcceptTemp;
      nReject += nRejectTemp;
    }
  }
}

/***************************************************************************
 * $RCSfile: DMCNonLocalUpdate.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

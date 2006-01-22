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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMCPbyPUpdate.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/CommCreate.h"

namespace qmcplusplus { 

  /// Constructor.
  DMCPbyPUpdate::DMCPbyPUpdate(MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& h,
      RandomGenerator_t& rg):
    W(w),Psi(psi),H(h), RandomGen(rg)
    { }
  
  /// destructor
  DMCPbyPUpdate::~DMCPbyPUpdate() { }

  void DMCPbyPUpdate::resetRun(BranchEngineType* brancher) {
    branchEngine=brancher;
    NumPtcl = W.getTotalNum();
    deltaR.resize(NumPtcl);
    G.resize(NumPtcl);
    dG.resize(NumPtcl);
    L.resize(NumPtcl);
    dL.resize(NumPtcl);

    Tau=brancher->Tau;
    m_oneover2tau = 1.0/(2.0*Tau);
    m_sqrttau = sqrt(Tau);
  }

  void DMCPbyPUpdate::resetBlock() {
    nAccept = 0; 
    nReject=0;
    nAllRejected = 0;
    nNodeCrossing=0;
  }

  /** advance all the walkers with killnode==yes
   * @param nat number of particles to move
   * 
   * When killnode==yes, any move resulting in node-crossing will lead to
   * the death of the walker.
   */
  void DMCPbyPUpdate::advanceKillNodeCrossing(WalkerIter_t it, WalkerIter_t it_end) {
    while(it != it_end) {

      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Buffer_t& w_buffer(thisWalker.DataSet);

      thisWalker.Weight = 1.0e0;
      thisWalker.Multiplicity=1.0e0;
      //save old local energy
      ValueType eold(thisWalker.Properties(LOCALENERGY));
      ValueType emixed(eold), enew(eold);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      //create a 3N-Dimensional Gaussian with variance=1
      //makeGaussRandom(deltaR);
      makeGaussRandomWithEngine(deltaR,Random);
      bool notcrossed(true);
      int nAcceptTemp(0);
      int nRejectTemp(0);
      int iat=0;
      RealType rr_proposed=0.0;
      RealType rr_accepted=0.0;
      while(notcrossed && iat<NumPtcl){

        PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        PosType newpos(W.makeMove(iat,dr));
        RealType ratio(Psi.ratio(W,iat,dG,dL));

        if(ratio < 0.0) {//node is crossed, stop here
          notcrossed = false;
          Psi.restore(iat);
          ++nNodeCrossing;
        } else {
          RealType rr=dot(dr,dr);
          rr_proposed+=rr;
          G = W.G+dG;
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
          ValueType scale=Tau;
          //ValueType vsq = dot(G[iat],G[iat]);
          //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
          dr = thisWalker.R[iat]-newpos-scale*G[iat]; 
          RealType logGb = -m_oneover2tau*dot(dr,dr);
          RealType prob = std::min(1.0,ratio*ratio*exp(logGb-logGf));
          if(RandomGen() < prob) { 
            ++nAcceptTemp;
            W.acceptMove(iat);
            Psi.update2(W,iat);
            W.G = G;
            W.L += dL;
            thisWalker.Drift = scale*G;
            rr_accepted+=rr;
          } else {
            ++nRejectTemp; 
            Psi.restore(iat);
          }
        } 
        ++iat;
      }

      if(notcrossed) {
        if(nAcceptTemp) {//need to overwrite the walker properties
          w_buffer.rewind();
          W.copyToBuffer(w_buffer);
          ValueType psi = Psi.evaluate(W,w_buffer);
          thisWalker.R = W.R;
          RealType enew= H.evaluate(W);
          thisWalker.resetProperty(log(abs(psi)),psi,enew);
          H.saveProperty(thisWalker.getPropertyBase());
          emixed = (eold+enew)*0.5;
          eold=enew;
        } else {
          thisWalker.Age++;
          ++nAllRejected;
        }

        RealType tau_eff=Tau;
        ValueType M = branchEngine->branchGF(tau_eff,emixed,0.0);
        if(thisWalker.Age > 3) M = std::min(0.5,M);
        else if(thisWalker.Age > 0) M = std::min(1.0,M);
        thisWalker.Weight = M; 
        thisWalker.Multiplicity=M + RandomGen();

        //accumulate the energy
        nAccept += nAcceptTemp;
        nReject += nRejectTemp;
      } else {//set the weight and multiplicity to zero
        thisWalker.willDie();
        nReject += NumPtcl;
      }
      branchEngine->accumulate(eold,thisWalker.Weight);
      ++it; 
    }
  }

  /** advance all the walkers with killnode==no
   * @param nat number of particles to move
   * 
   * When killnode==no, any move resulting in node-crossing is treated
   * as a normal rejection.
   */
  void DMCPbyPUpdate::advanceRejectNodeCrossing(WalkerIter_t it, WalkerIter_t it_end) {
    while(it != it_end) {

      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Buffer_t& w_buffer(thisWalker.DataSet);

      thisWalker.Weight = 1.0e0;
      thisWalker.Multiplicity=1.0e0;
      //save old local energy
      ValueType eold(thisWalker.Properties(LOCALENERGY));
      ValueType emixed(eold), enew(eold);

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

        if(ratio < 0.0) {//node is crossed reject the move
          ++nRejectTemp;
          ++nNodeCrossing;
          Psi.restore(iat);
        } else {
          G = W.G+dG;
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
          RealType scale=Tau;
          //ValueType vsq = Dot(G,G);
          //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
          dr = thisWalker.R[iat]-newpos-scale*G[iat]; 
          RealType logGb = -m_oneover2tau*dot(dr,dr);
          RealType prob = std::min(1.0,ratio*ratio*exp(logGb-logGf));
          if(RandomGen() < prob) { 
            ++nAcceptTemp;
            W.acceptMove(iat);
            Psi.update2(W,iat);
            W.G = G;
            W.L += dL;
            thisWalker.Drift = scale*G;
            rr_accepted+=rr;
          } else {
            ++nRejectTemp; 
            Psi.restore(iat);
          }
        } 

        ++iat;
      }

      if(nAcceptTemp>0) {//need to overwrite the walker properties
        thisWalker.R = W.R;
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        ValueType psi = Psi.evaluate(W,w_buffer);
        enew= H.evaluate(W);
        thisWalker.resetProperty(log(abs(psi)),psi,enew);
        H.saveProperty(thisWalker.getPropertyBase());
        emixed = (eold+enew)*0.5e0;
      } else {
        thisWalker.Age++;
        ++nAllRejected;
        rr_accepted=0.0;
      }

      ValueType M = branchEngine->branchGF(Tau*rr_accepted/rr_proposed,emixed,0.0);
      if(thisWalker.Age > 1) M = std::min(0.5,M);
      else if(thisWalker.Age > 0) M = std::min(1.0,M);
      thisWalker.Weight = M; 
      thisWalker.Multiplicity=M + RandomGen();
      branchEngine->accumulate(emixed,thisWalker.Weight);//accumulate the energy
      nAccept += nAcceptTemp;
      nReject += nRejectTemp;
      ++it;
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

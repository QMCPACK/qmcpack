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
#include "Message/CommCreate.h"

namespace qmcplusplus { 

  /// Constructor.
  DMCUpdatePbyPWithRejection::DMCUpdatePbyPWithRejection(ParticleSet& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      RandomGenerator_t& rg): QMCUpdateBase(w,psi,h,rg)
    { }
  
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
    int item=0;
    while(it != it_end) {
      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      //thisWalker.Weight = 1.0e0;
      //thisWalker.Multiplicity=1.0e0;
      //save old local energy
      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType emixed(eold), enew(eold);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(deltaR,RandomGen);
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
          RealType scale=Tau;
          //ValueType vsq = Dot(G,G);
          //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
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
      }

      if(nAcceptTemp>0) {//need to overwrite the walker properties
        thisWalker.R = W.R;
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        RealType psi = Psi.evaluate(W,w_buffer);
        enew= H.evaluate(W);
        thisWalker.resetProperty(std::log(abs(psi)),psi,enew);
        H.saveProperty(thisWalker.getPropertyBase());
        emixed = (eold+enew)*0.5e0;
      } else {
        thisWalker.Age++;
        ++nAllRejected;
        rr_accepted=0.0;
      }

      //branchEngine->setWeight(thisWalker,Tau*rr_accepted/rr_proposed,emixed,RandomGen());
      thisWalker.Weight *= branchEngine->branchGF(Tau*rr_accepted/rr_proposed,emixed,0.0);
      //if(MaxAge) {
      //  RealType M=thisWalker.Weight;
      //  if(thisWalker.Age > MaxAge) M = std::min(0.5,M);
      //  else if(thisWalker.Age > 0) M = std::min(1.0,M);
      //  thisWalker.Multiplicity = M + RandomGen();
      //  branchEngine->accumulate(eold,M);
      //} else {
      //  branchEngine->accumulate(eold,1);
      //}
      //branchEngine->accumulate(eold,1);
      
      nAccept += nAcceptTemp;
      nReject += nRejectTemp;
      ++it;
    }
  }

  /// Constructor.
  DMCUpdatePbyPWithKill::DMCUpdatePbyPWithKill(ParticleSet& w, TrialWaveFunction& psi, QMCHamiltonian& h,
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

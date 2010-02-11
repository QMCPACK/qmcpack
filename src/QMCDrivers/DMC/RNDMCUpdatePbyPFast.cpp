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


namespace qmcplusplus { 

  /// Constructor.
  RNDMCUpdatePbyPWithRejectionFast::RNDMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg): 
    QMCUpdateBase(w,psi,h,rg)
    { 
      myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::advance")); //timer for the walker loop
      myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::movePbyP")); //timer for MC, ratio etc
      myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::updateMBO")); //timer for measurements
      myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::energy")); //timer for measurements
      TimerManager.addTimer(myTimers[0]);
      TimerManager.addTimer(myTimers[1]);
      TimerManager.addTimer(myTimers[2]);
      TimerManager.addTimer(myTimers[3]);
    }
  
  /// destructor
  RNDMCUpdatePbyPWithRejectionFast::~RNDMCUpdatePbyPWithRejectionFast() { }

  /** advance all the walkers with killnode==no
   * @param nat number of particles to move
   * 
   * When killnode==no, any move resulting in node-crossing is treated
   * as a normal rejection.
   */
  void RNDMCUpdatePbyPWithRejectionFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end
      , bool measure) 
  {
    myTimers[0]->start();
    for(;it != it_end;++it) 
    {
      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      W.loadWalker(thisWalker,true);
      //W.R = thisWalker.R;
      //w_buffer.rewind();
      //W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);
      int oldRNsign = nNodeCrossing;

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
        //get the displacement
        GradType grad_iat=Psi.evalGrad(W,iat);
        PosType dr;
        //RealType sc=getDriftScale(m_tauovermass,grad_iat);
        //PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_iat));
        getScaledDrift(m_tauovermass, grad_iat, dr);
        dr += m_sqrttau * deltaR[iat];

        //RealType rr=dot(dr,dr);
        RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
        rr_proposed+=rr;

        if(rr>m_r2max)
        {
          ++nRejectTemp; continue;
        }

        //PosType newpos(W.makeMove(iat,dr));
        if(!W.makeMoveAndCheck(iat,dr)) continue;
        PosType newpos(W.R[iat]);
        RealType ratio = Psi.ratioGrad(W,iat,grad_iat);
        bool valid_move=false;

        //node is crossed reject the move
        //if(Psi.getPhase() > numeric_limits<RealType>::epsilon()) 
        //if(branchEngine->phaseChanged(Psi.getPhase(),thisWalker.Properties(SIGN))) 
        if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
        {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.rejectMove(iat); Psi.rejectMove(iat);
        } 
        else 
        {
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);

          //Use the force of the particle iat
          //RealType scale=getDriftScale(m_tauovermass,grad_iat);
          //dr = thisWalker.R[iat]-newpos-scale*real(grad_iat);
          getScaledDrift(m_tauovermass, grad_iat, dr);
          dr = thisWalker.R[iat] - newpos - dr;

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
            
          } 
          else 
          {
            ++nRejectTemp; 
            W.rejectMove(iat); Psi.rejectMove(iat);
          }
        } 
      }
      myTimers[1]->stop();

      //RealType nodecorr_old=thisWalker.Properties(DRIFTSCALE);
      //RealType nodecorr=nodecorr_old; 
      bool advanced=true;

      if(nAcceptTemp>0) 
      {//need to overwrite the walker properties
        myTimers[2]->start();
        thisWalker.Age=0;
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
      {//all moves are rejected: does not happen normally with reasonable wavefunctions
        advanced=false;
        thisWalker.Age++;
        thisWalker.Properties(R2ACCEPTED)=0.0;
        H.rejectedMove(W,thisWalker); 
        ++nAllRejected;
        enew=eold;//copy back old energy
        gf_acc=1.0;
      }
      
      
//       RealType sigma_last = thisWalker.ReleasedNodeWeight/std::abs(thisWalker.ReleasedNodeWeight);
      thisWalker.Weight *= branchEngine->branchWeightBare(enew,eold);
      ValueType altR = Psi.alternateRatio(W);
      thisWalker.ReleasedNodeWeight = altR/std::abs(altR);
      RealType KE_ferm_n = -0.5*(Sum(W.L)+Dot(W.G,W.G)); 
      thisWalker.resetReleasedNodeProperty(enew,thisWalker.Properties(LOCALENERGY)+KE_ferm_n);
      nAccept += nAcceptTemp;
      nReject += nRejectTemp;

    }
    myTimers[0]->stop();
  }

}

/***************************************************************************
 * $RCSfile: DMCUpdatePbyP.cpp,v $   $Author: mcminis2 $
 * $Revision: 3697 $   $Date: 2009-03-24 19:30:46 -0400 (Tue, 24 Mar 2009) $
 * $Id: DMCUpdatePbyP.cpp 3697 2009-03-24 23:30:46Z mcminis2 $ 
 ***************************************************************************/

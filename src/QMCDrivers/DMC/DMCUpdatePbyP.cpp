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
//#define TEST_INNERBRANCH

namespace qmcplusplus { 

  /// Constructor.
  DMCUpdatePbyPWithRejection::DMCUpdatePbyPWithRejection(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      RandomGenerator_t& rg): QMCUpdateBase(w,psi,h,rg)
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

#if defined(TEST_INNERBRANCH)
    for(;it != it_end;++it) 
    {
      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      for(int step=0; step<5; ++step)
      {
        //create a 3N-Dimensional Gaussian with variance=1
        makeGaussRandomWithEngine(deltaR,RandomGen);
        int nAcceptTemp(0);
        int nRejectTemp(0);
        RealType eold(thisWalker.Properties(LOCALENERGY));
        RealType enew(eold);
        RealType rr_proposed=0.0;
        RealType rr_accepted=0.0;
        for(int iat=0; iat<NumPtcl; ++iat) 
        {
          //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
          RealType sc=getDriftScale(Tau,W.G[iat]);
          PosType dr(m_sqrttau*deltaR[iat]+sc*W.G[iat]);

          //RealType rr=dot(dr,dr);
          RealType rr=Tau*dot(deltaR[iat],deltaR[iat]);
          rr_proposed+=rr;

          if(rr>m_r2max)//reject a big move
          {
            ++nRejectTemp; continue;
          }

          PosType newpos(W.makeMove(iat,dr));
          RealType ratio=Psi.ratio(W,iat,dG,dL);

          ///node is crossed reject the move
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

            //Scale is set by the total quantum force
            //RealType scale=getDriftScale(Tau,G);
            //dr = thisWalker.R[iat]-newpos-scale*real(G[iat]); 
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

              //Checking pbyp drift
              //assignDrift(scale,G,thisWalker.Drift);

              rr_accepted+=rr;
            } else {
              ++nRejectTemp; 
              W.rejectMove(iat); Psi.rejectMove(iat);
            }
          } 
        }

        if(nAcceptTemp>0) 
        {//need to overwrite the walker properties
          thisWalker.R = W.R;
          w_buffer.rewind();
          W.copyToBuffer(w_buffer);
          RealType psi = Psi.evaluate(W,w_buffer);
          enew= H.evaluate(W);
          thisWalker.resetProperty(std::log(abs(psi)),psi,enew,rr_accepted,rr_proposed,1.0);
          H.saveProperty(thisWalker.getPropertyBase());

          //update the drift: safe operator for QMC_COMPLEX=1
          PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
          //thisWalker.Drift=W.G;
          //setScaledDriftPbyP(Tau,W.G,thisWalker.Drift);
        } else {
          thisWalker.Age++;
          thisWalker.Properties(R2ACCEPTED)=0.0;
          enew=eold;//copy back old energy
          ++nAllRejected;
        }

        thisWalker.Weight *= branchEngine->branchWeight(eold,enew);
        eold=enew;
        nAccept += nAcceptTemp;
        nReject += nRejectTemp;
      }
    }
#else
    myTimers[0]->start();
    for(;it != it_end;++it) 
    {
      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
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

      myTimers[1]->start();
      for(int iat=0; iat<NumPtcl; ++iat) 
      {

        //get the displacement
        //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        RealType sc=getDriftScale(Tau,W.G[iat]);
        PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));

        //RealType rr=dot(dr,dr);
        RealType rr=Tau*dot(deltaR[iat],deltaR[iat]);
        rr_proposed+=rr;

        if(rr>m_r2max)
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

          //Scale is set by the total quantum force
          //RealType scale=getDriftScale(Tau,G);
          //dr = thisWalker.R[iat]-newpos-scale*real(G[iat]); 
          
          //Use the force of the particle iat
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
      }
      myTimers[1]->stop();
      
      if(nAcceptTemp>0) 
      {//need to overwrite the walker properties
        myTimers[2]->start();
        thisWalker.R = W.R;

        //copy the new Gradient to drift
        PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
        //setScaledDriftPbyP(Tau,W.G,(*it)->Drift);

        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        RealType psi = Psi.evaluate(W,w_buffer);
        myTimers[2]->stop();

        myTimers[3]->start();
        enew= H.evaluate(W);
        myTimers[3]->stop();

        thisWalker.resetProperty(std::log(abs(psi)),psi,enew,rr_accepted,rr_proposed,1.0);
        H.saveProperty(thisWalker.getPropertyBase());
      } 
      else 
      {
        thisWalker.Age++;
        thisWalker.Properties(R2ACCEPTED)=0.0;
        ++nAllRejected;
        enew=eold;//copy back old energy
      }

      thisWalker.Weight *= branchEngine->branchWeight(eold,enew);

      nAccept += nAcceptTemp;
      nReject += nRejectTemp;

    }
    myTimers[0]->stop();
#endif
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

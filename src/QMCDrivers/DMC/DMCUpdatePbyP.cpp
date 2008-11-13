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

//define macros to print out runtime data
#if defined(PRINT_DEBUG)
#define DMC_TRACE_START(NOW) NOW
#define DMC_TRACE_STOP(WID,PID,MVD,ELAPSED) \
  OhmmsInfo::Debug->getStream() << setw(16) << WID \
  << setw(5) << PID << setw(4) << MVD << setw(15) << ELAPSED << std::endl
#else
#define DMC_TRACE_START(NOW) 
#define DMC_TRACE_STOP(WID,MID,PID,ELAPSED) 
#endif

namespace qmcplusplus { 

  /// Constructor.
  DMCUpdatePbyPWithRejectionFast::DMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg): 
    QMCUpdateBase(w,psi,h,rg)
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
  DMCUpdatePbyPWithRejectionFast::~DMCUpdatePbyPWithRejectionFast() { }

  /** advance all the walkers with killnode==no
   * @param nat number of particles to move
   * 
   * When killnode==no, any move resulting in node-crossing is treated
   * as a normal rejection.
   */
  void DMCUpdatePbyPWithRejectionFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
      bool measure) 
  {

    Timer localTimer;
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
      RealType gf_acc=1.0;

      myTimers[1]->start();
      for(int iat=0; iat<NumPtcl; ++iat) 
      {

        DMC_TRACE_START(localTimer.restart());

        //get the displacement
        //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        GradType grad_iat=Psi.evalGrad(W,iat);
        RealType sc=getDriftScale(m_tauovermass,grad_iat);
        PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_iat));

        //RealType rr=dot(dr,dr);
        RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
        rr_proposed+=rr;

        if(rr>m_r2max)
        {
          ++nRejectTemp; continue;
        }

        PosType newpos(W.makeMove(iat,dr));
        RealType ratio = Psi.ratioGrad(W,iat,grad_iat);
        bool valid_move=false;

        //node is crossed reject the move
        if(Psi.getPhase() > numeric_limits<RealType>::epsilon()) 
        {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.rejectMove(iat); Psi.rejectMove(iat);
        } 
        else 
        {
          RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);

          //Use the force of the particle iat
          RealType scale=getDriftScale(m_tauovermass,grad_iat);
          dr = thisWalker.R[iat]-newpos-scale*real(grad_iat); 

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

        DMC_TRACE_STOP(thisWalker.ID,iat,valid_move,localTimer.elapsed());
      }
      myTimers[1]->stop();

      DMC_TRACE_START(localTimer.restart());
      
      //RealType nodecorr_old=thisWalker.Properties(DRIFTSCALE);
      //RealType nodecorr=nodecorr_old;
      bool advanced=true;

      if(nAcceptTemp>0) 
      {//need to overwrite the walker properties
        myTimers[2]->start();
        thisWalker.Age=0;
        thisWalker.R = W.R;

        w_buffer.rewind();
        W.updateBuffer(w_buffer);
        RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
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
        thisWalker.rejectedMove();
        ++nAllRejected;
        enew=eold;//copy back old energy
        gf_acc=1.0;
      }

      //2008-06-26: select any
      //bare green function by setting nodecorr=nodecorr_old=1.0
      thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
      DMC_TRACE_STOP(thisWalker.ID,NumPtcl,advanced,localTimer.elapsed());

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

//#if defined(TEST_INNERBRANCH)
//    for(;it != it_end;++it) 
//    {
//      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
//      Walker_t& thisWalker(**it);
//      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//
//      W.R = thisWalker.R;
//      w_buffer.rewind();
//      W.copyFromBuffer(w_buffer);
//      Psi.copyFromBuffer(W,w_buffer);
//
//      for(int step=0; step<5; ++step)
//      {
//        //create a 3N-Dimensional Gaussian with variance=1
//        makeGaussRandomWithEngine(deltaR,RandomGen);
//        int nAcceptTemp(0);
//        int nRejectTemp(0);
//        RealType eold(thisWalker.Properties(LOCALENERGY));
//        RealType enew(eold);
//        RealType rr_proposed=0.0;
//        RealType rr_accepted=0.0;
//        for(int iat=0; iat<NumPtcl; ++iat) 
//        {
//          //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
//          RealType sc=getDriftScale(Tau,W.G[iat]);
//          PosType dr(m_sqrttau*deltaR[iat]+sc*W.G[iat]);
//
//          //RealType rr=dot(dr,dr);
//          RealType rr=Tau*dot(deltaR[iat],deltaR[iat]);
//          rr_proposed+=rr;
//
//          if(rr>m_r2max)//reject a big move
//          {
//            ++nRejectTemp; continue;
//          }
//
//          PosType newpos(W.makeMove(iat,dr));
//          RealType ratio=Psi.ratio(W,iat,dG,dL);
//
//          ///node is crossed reject the move
//          if(Psi.getPhase() > numeric_limits<RealType>::epsilon()) 
//          {
//            ++nRejectTemp;
//            ++nNodeCrossing;
//            W.rejectMove(iat); Psi.rejectMove(iat);
//          } 
//          else 
//          {
//            G = W.G+dG;
//            RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
//
//            //Scale is set by the total quantum force
//            //RealType scale=getDriftScale(Tau,G);
//            //dr = thisWalker.R[iat]-newpos-scale*real(G[iat]); 
//            RealType scale=getDriftScale(Tau,G[iat]);
//            dr = thisWalker.R[iat]-newpos-scale*real(G[iat]); 
//
//            RealType logGb = -m_oneover2tau*dot(dr,dr);
//            RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
//            if(RandomGen() < prob) 
//            { 
//              ++nAcceptTemp;
//              W.acceptMove(iat);
//              Psi.acceptMove(W,iat);
//              W.G = G;
//              W.L += dL;
//
//              //Checking pbyp drift
//              //assignDrift(scale,G,thisWalker.Drift);
//
//              rr_accepted+=rr;
//            } else {
//              ++nRejectTemp; 
//              W.rejectMove(iat); Psi.rejectMove(iat);
//            }
//          } 
//        }
//
//        if(nAcceptTemp>0) 
//        {//need to overwrite the walker properties
//          thisWalker.R = W.R;
//          w_buffer.rewind();
//          W.copyToBuffer(w_buffer);
//          RealType psi = Psi.evaluate(W,w_buffer);
//          enew= H.evaluate(W);
//          thisWalker.resetProperty(std::log(abs(psi)),psi,enew,rr_accepted,rr_proposed,1.0);
//          H.saveProperty(thisWalker.getPropertyBase());
//
//          //update the drift: safe operator for QMC_COMPLEX=1
//          //2008-08-26 THIS IS NOT DOING ANYTHING
//          //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
//        } else {
//          thisWalker.Age++;
//          thisWalker.Properties(R2ACCEPTED)=0.0;
//          enew=eold;//copy back old energy
//          ++nAllRejected;
//        }
//
//        thisWalker.Weight *= branchEngine->branchWeight(eold,enew);
//        eold=enew;
//        nAccept += nAcceptTemp;
//        nReject += nRejectTemp;
//      }
//    }
//#else
    Timer localTimer;

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
      RealType gf_acc=1.0;

      myTimers[1]->start();
      for(int iat=0; iat<NumPtcl; ++iat) 
      {

        DMC_TRACE_START(localTimer.restart());

        //get the displacement
        //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
        RealType sc=getDriftScale(m_tauovermass,W.G[iat]);
        PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));

        //RealType rr=dot(dr,dr);
        RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
        rr_proposed+=rr;

        if(rr>m_r2max)
        {
//           cout<<" PROPOSED move was too big!!"<<endl;
          ++nRejectTemp; continue;
        }

        PosType newpos(W.makeMove(iat,dr));
        RealType ratio=Psi.ratio(W,iat,dG,dL);
        bool valid_move=false;

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
          RealType scale=getDriftScale(m_tauovermass,G[iat]);
          dr = thisWalker.R[iat]-newpos-scale*real(G[iat]); 

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
            W.rejectMove(iat); Psi.rejectMove(iat);
          }
        } 

        DMC_TRACE_STOP(thisWalker.ID,iat,valid_move,localTimer.elapsed());
      }
      myTimers[1]->stop();

      DMC_TRACE_START(localTimer.restart());
      
      RealType nodecorr_old=thisWalker.Properties(DRIFTSCALE);
      RealType nodecorr=nodecorr_old;
      bool advanced=true;

      if(nAcceptTemp>0) 
      {//need to overwrite the walker properties
        myTimers[2]->start();
        thisWalker.Age=0;
        thisWalker.R = W.R;

        nodecorr=getNodeCorrection(W.G,thisWalker.Drift);

        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        //RealType psi = Psi.evaluate(W,w_buffer);
        RealType logpsi = Psi.evaluateLog(W,w_buffer);
        myTimers[2]->stop();

        myTimers[3]->start();
        enew= H.evaluate(W);
        myTimers[3]->stop();

        //thisWalker.resetProperty(std::log(abs(psi)),psi,enew,rr_accepted,rr_proposed,nodecorr);
        thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr );
        H.auxHevaluate(W,thisWalker);
        H.saveProperty(thisWalker.getPropertyBase());
      } 
      else 
      {//all moves are rejected: does not happen normally with reasonable wavefunctions
        advanced=false;
        thisWalker.rejectedMove();
        thisWalker.Age++;
        thisWalker.Properties(R2ACCEPTED)=0.0;
        ++nAllRejected;
        enew=eold;//copy back old energy
        gf_acc=1.0;
      }

      //2008-06-26: select any
      //bare green function by setting nodecorr=nodecorr_old=1.0
      thisWalker.Weight *= branchEngine->branchWeight(enew,eold);

      DMC_TRACE_STOP(thisWalker.ID,NumPtcl,advanced,localTimer.elapsed());

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
//#endif
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

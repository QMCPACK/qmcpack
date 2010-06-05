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
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

/** add timers for VMC PbyP updates
 * @param timers container for the timers
 */
void add_vmc_timers(vector<NewTimer*>& timers)
{
  timers.push_back(new NewTimer("VMCUpdatePbyP::advance")); //timer for the walker loop
  timers.push_back(new NewTimer("VMCUpdatePbyP::movePbyP")); //timer for MC, ratio etc
  timers.push_back(new NewTimer("VMCUpdatePbyP::updateMBO")); //timer for measurements
  timers.push_back(new NewTimer("VMCUpdatePbyP::energy")); //timer for measurements
  for (int i=0; i<timers.size(); ++i) TimerManager.addTimer(timers[i]);
}

VMCUpdatePbyP::VMCUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                             QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg)
{
  add_vmc_timers(myTimers);
}

VMCUpdatePbyP::~VMCUpdatePbyP()
{
}

void VMCUpdatePbyP::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{

  myTimers[0]->start();
  for (; it != it_end; ++it)
    {

      Walker_t& thisWalker(**it);
      W.loadWalker(thisWalker,true);

      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      Psi.copyFromBuffer(W,w_buffer);

      myTimers[1]->start();
      for (int iter=0; iter<nSubSteps; ++iter)
        {
          makeGaussRandomWithEngine(deltaR,RandomGen);
          bool stucked=true;
          for (int iat=0; iat<W.getTotalNum(); ++iat)
            {
              PosType dr = m_sqrttau*deltaR[iat];
              //replace makeMove by makeMoveAndCheck
              //PosType newpos = W.makeMove(iat,dr);
              if (W.makeMoveAndCheck(iat,dr))
                {
                  RealType ratio = Psi.ratio(W,iat);
                  RealType prob = ratio*ratio;
                  //RealType prob = std::min(1.0e0,ratio*ratio);
                  if (RandomGen() < prob)
                    {
                      stucked=false;
                      ++nAccept;
                      W.acceptMove(iat);
                      Psi.acceptMove(W,iat);
                    }
                  else
                    {
                      ++nReject;
                      W.rejectMove(iat);
                      Psi.rejectMove(iat);
                    }
                }
              else //reject illegal moves
                ++nReject;
            }

          if (stucked)
            {
              ++nAllRejected;
              H.rejectedMove(W,thisWalker);
            }
        }
      myTimers[1]->stop();

      myTimers[2]->start();
      //thisWalker.R = W.R;
      //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
      //w_buffer.rewind();
      //W.updateBuffer(w_buffer);

      RealType logpsi = Psi.updateBuffer(W,w_buffer,true);
      W.saveWalker(thisWalker);

      //W.copyToBuffer(w_buffer);
      //RealType logpsi = Psi.evaluate(W,w_buffer);
      myTimers[2]->stop();

      myTimers[3]->start();
      RealType eloc=H.evaluate(W);
      myTimers[3]->stop();

      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());

    }
  myTimers[0]->stop();
}

/// Constructor.
VMCUpdatePbyPWithDrift::VMCUpdatePbyPWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg)
{
  add_vmc_timers(myTimers);
}

VMCUpdatePbyPWithDrift::~VMCUpdatePbyPWithDrift()
{
}

void VMCUpdatePbyPWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  myTimers[0]->start();
  for (; it != it_end; ++it)
    {
      Walker_t& thisWalker(**it);

      W.loadWalker(thisWalker,true);

      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      Psi.copyFromBuffer(W,thisWalker.DataSet);
      //create a 3N-Dimensional Gaussian with variance=1
      bool moved = false;
      myTimers[1]->start();
      for (int iter=0; iter<nSubSteps; ++iter)
        {

          makeGaussRandomWithEngine(deltaR,RandomGen);

          moved = false;

          for (int iat=0; iat<W.getTotalNum(); ++iat)
            {
              PosType dr;
              ///dr = m_sqrttau*deltaR[iat]+thisWalker.Drift[iat];
              //RealType sc=getDriftScale(m_tauovermass,W.G[iat]);
              //PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));
              getScaledDrift(m_tauovermass,W.G[iat],dr);
              dr += m_sqrttau*deltaR[iat];

              //reject illegal moves
              if (!W.makeMoveAndCheck(iat,dr))
                {
                  ++nReject;
                  continue;
                }
              //PosType newpos=W.R[iat];
              //PosType newpos = W.makeMove(iat,dr);

              RealType ratio = Psi.ratio(W,iat,dG,dL);
              RealType prob = ratio*ratio;

              //zero is always rejected
              if (prob<numeric_limits<RealType>::epsilon())
                {
                  ++nReject;
                  W.rejectMove(iat);
                  Psi.rejectMove(iat);
                  continue;
                }

              G = W.G+dG;

              //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
              //dr = (*it)->R[iat]-newpos-Tau*G[iat];
              //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
              RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);

              //RealType scale=getDriftScale(m_tauovermass,G[iat]);
              //dr = thisWalker.R[iat]-W.R[iat]-scale*real(G[iat]);
              getScaledDrift(m_tauovermass,G[iat],dr);
              dr = thisWalker.R[iat]-W.R[iat]-dr;

              RealType logGb = -m_oneover2tau*dot(dr,dr);

              //RealType prob = std::min(1.0e0,ratio*ratio*std::exp(logGb-logGf));
              if (RandomGen() < prob*std::exp(logGb-logGf))
                {
                  moved = true;
                  ++nAccept;
                  W.acceptMove(iat);
                  Psi.acceptMove(W,iat);
                  W.G = G;
                  W.L += dL;

                  //do not need to update Drift
                  //assignDrift(scale,G,thisWalker.Drift);

                }
              else
                {
                  ++nReject;
                  W.rejectMove(iat);
                  Psi.rejectMove(iat);
                }
            }
        }
      myTimers[1]->stop();

      if (moved)
        {
          myTimers[2]->start();
          //thisWalker.R = W.R;
          //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
          //w_buffer.rewind();
          //W.copyToBuffer(w_buffer);

          RealType logpsi = Psi.evaluateLog(W,w_buffer);
          W.saveWalker(thisWalker);

          myTimers[2]->stop();

          myTimers[3]->start();
          RealType eloc=H.evaluate(W);
          myTimers[3]->stop();
          //thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
          thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
          H.auxHevaluate(W,thisWalker);
          H.saveProperty(thisWalker.getPropertyBase());
        }
      else
        {
          ++nAllRejected;
          H.rejectedMove(W,thisWalker);
        }
    }
  myTimers[0]->stop();
}

/// Constructor.
VMCUpdatePbyPWithDriftFast::VMCUpdatePbyPWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg) :
    QMCUpdateBase(w,psi,h,rg)
{
  add_vmc_timers(myTimers);
}

VMCUpdatePbyPWithDriftFast::~VMCUpdatePbyPWithDriftFast()
{
}

void VMCUpdatePbyPWithDriftFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  myTimers[0]->start();
  for (; it != it_end; ++it)
    {
      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      W.loadWalker(thisWalker,true);
      //W.R = thisWalker.R;
      //w_buffer.rewind();
      //W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);


      myTimers[1]->start();
      bool moved = false;
      for (int iter=0; iter<nSubSteps; ++iter)
        {
          //create a 3N-Dimensional Gaussian with variance=1
          makeGaussRandomWithEngine(deltaR,RandomGen);
          moved = false;
          for (int iat=0; iat<W.getTotalNum(); ++iat)
            {

              GradType grad_now=Psi.evalGrad(W,iat), grad_new;
              PosType dr;
//    //RealType sc=getDriftScale(m_tauovermass,grad_now);
//    //PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_now));
              getScaledDrift(m_tauovermass,grad_now,dr);
              dr += m_sqrttau*deltaR[iat];
//    PosType dr(m_sqrttau*deltaR[iat]+m_tauovermass*grad_now);
              if (!W.makeMoveAndCheck(iat,dr))
                {
                  ++nReject;
                  continue;
                }

              //PosType newpos = W.makeMove(iat,dr);
              RealType ratio = Psi.ratioGrad(W,iat,grad_new);
              RealType prob = ratio*ratio;

              //zero is always rejected
              if (prob<numeric_limits<RealType>::epsilon())
                {
                  ++nReject;
                  W.rejectMove(iat);
                  Psi.rejectMove(iat);
                  continue;
                }

              //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
              //dr = (*it)->R[iat]-newpos-Tau*G[iat];
              //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
              RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);

//    //sc=getDriftScale(m_tauovermass,grad_new);
//    //dr = thisWalker.R[iat]-W.R[iat]-sc*real(grad_new);
              getScaledDrift(m_tauovermass,grad_new,dr);
              dr = thisWalker.R[iat]-W.R[iat]-dr;
//    dr = thisWalker.R[iat]-W.R[iat]-m_tauovermass*grad_new;
              RealType logGb = -m_oneover2tau*dot(dr,dr);

              //RealType prob = std::min(1.0e0,ratio*ratio*std::exp(logGb-logGf));
              if (RandomGen() < prob*std::exp(logGb-logGf))
                {
                  moved = true;
                  ++nAccept;
                  W.acceptMove(iat);
                  Psi.acceptMove(W,iat);
                }
              else
                {
                  ++nReject;
                  W.rejectMove(iat);
                  Psi.rejectMove(iat);
                }
            }
          //for subSteps must update thiswalker
          thisWalker.R=W.R;
          thisWalker.G=W.G;
          thisWalker.L=W.L;
          myTimers[1]->stop();
        }

      if (moved)
        {
          myTimers[2]->start();
          //thisWalker.R = W.R;
          //w_buffer.rewind();
          //W.updateBuffer(w_buffer);
          RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
          W.saveWalker(thisWalker);
          myTimers[2]->stop();

          myTimers[3]->start();
          RealType eloc=H.evaluate(W);
          myTimers[3]->stop();
          //thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
          thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
          H.auxHevaluate(W,thisWalker);
          H.saveProperty(thisWalker.getPropertyBase());
        }
      else
        {
          ++nAllRejected;
          H.rejectedMove(W,thisWalker);
        }

    }
  myTimers[0]->stop();
}


VMCUpdatePbyPSamplePsi::VMCUpdatePbyPSamplePsi(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg)
{
  add_vmc_timers(myTimers);
}

VMCUpdatePbyPSamplePsi::~VMCUpdatePbyPSamplePsi()
{
}

void VMCUpdatePbyPSamplePsi::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{

  myTimers[0]->start();
  for (; it != it_end; ++it)
    {

      Walker_t& thisWalker(**it);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

      W.loadWalker(thisWalker,true);
      //W.R = thisWalker.R;
      //w_buffer.rewind();
      //W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      myTimers[1]->start();
      for (int iter=0; iter<nSubSteps; ++iter)
        {
          makeGaussRandomWithEngine(deltaR,RandomGen);
          bool stucked=true;
          for (int iat=0; iat<W.getTotalNum(); ++iat)
            {
              PosType dr = m_sqrttau*deltaR[iat];
              //ignore illegal moves
              if (!W.makeMoveAndCheck(iat,dr)) continue;
              //PosType newpos = W.makeMove(iat,dr);

              RealType ratio = Psi.ratio(W,iat);
              RealType prob = ratio;
              //RealType prob = std::min(1.0e0,ratio*ratio);
              if (RandomGen() < prob)
                {
                  stucked=false;
                  ++nAccept;
                  W.acceptMove(iat);
                  Psi.acceptMove(W,iat);
                }
              else
                {
                  ++nReject;
                  W.rejectMove(iat);
                  Psi.rejectMove(iat);
                }
            }
          if (stucked)
            {
              ++nAllRejected;
              H.rejectedMove(W,thisWalker);
            }
        }
      myTimers[1]->stop();

      myTimers[2]->start();
      //thisWalker.R = W.R;
      //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
      //w_buffer.rewind();
      //W.updateBuffer(w_buffer);
      RealType logpsi = Psi.updateBuffer(W,w_buffer,true);
      W.saveWalker(thisWalker);

      myTimers[2]->stop();

      myTimers[3]->start();
      RealType eloc=H.evaluate(W);
      myTimers[3]->stop();

      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.auxHevaluate(W,thisWalker);
      thisWalker.Weight=std::exp(logpsi);
      H.saveProperty(thisWalker.getPropertyBase());

    }
  myTimers[0]->stop();
}



}

/***************************************************************************
 * $RCSfile: VMCUpdatePbyP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdatePbyP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/

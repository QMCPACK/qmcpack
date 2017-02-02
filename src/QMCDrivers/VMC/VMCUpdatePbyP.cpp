//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

namespace qmcplusplus
{

/** add timers for VMC PbyP updates
 * @param timers container for the timers
 */
void add_vmc_timers(std::vector<NewTimer*>& timers)
{
  timers.push_back(new NewTimer("VMCUpdatePbyP::advance",timer_level_medium)); //timer for the walker loop
  timers.push_back(new NewTimer("VMCUpdatePbyP::movePbyP",timer_level_medium)); //timer for MC, ratio etc
  timers.push_back(new NewTimer("VMCUpdatePbyP::updateMBO",timer_level_medium)); //timer for measurements
  timers.push_back(new NewTimer("VMCUpdatePbyP::energy",timer_level_medium)); //timer for measurements
  for (int i=0; i<timers.size(); ++i)
    TimerManager.addTimer(timers[i]);
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

void VMCUpdatePbyP::advanceWalker(Walker_t& thisWalker)
{
  W.loadWalker(thisWalker,true);
  Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
  Psi.copyFromBuffer(W,w_buffer);
  myTimers[1]->start();
  for (int iter=0; iter<nSubSteps; ++iter)
  {
    makeGaussRandomWithEngine(deltaR,RandomGen);
    bool stucked=true;
    for(int ig=0; ig<W.groups(); ++ig) //loop over species
    {
      RealType sqrttau = std::sqrt(Tau*MassInvS[ig]);
      for (int iat=W.first(ig); iat<W.last(ig); ++iat)
      {
        mPosType dr = sqrttau*deltaR[iat];
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
      } //iat
    }//ig for the species
    if (stucked)
    {
      ++nAllRejected;
      //H.rejectedMove(W,thisWalker);
    }
    thisWalker.R=W.R;
  }
  myTimers[1]->stop();
  myTimers[2]->start();
  //thisWalker.R = W.R;
  //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
  //w_buffer.rewind();
  //W.updateBuffer(w_buffer);
  RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
  W.saveWalker(thisWalker);
  //W.copyToBuffer(w_buffer);
  //RealType logpsi = Psi.evaluate(W,w_buffer);
  myTimers[2]->stop();
  myTimers[3]->start();
  EstimatorRealType eloc=H.evaluate(W);
  myTimers[3]->stop();
  thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
  H.auxHevaluate(W,thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());
}

void VMCUpdatePbyP::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  myTimers[0]->start();
  for (; it != it_end; ++it)
  {
    advanceWalker(**it);
  }
  myTimers[0]->stop();
}

VMCUpdatePbyP::RealType VMCUpdatePbyP::advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios)
{
  Walker_t::Buffer_t& w_buffer(w1.DataSet);
  W.loadWalker(w1,true);
  Psi.copyFromBuffer(W,w_buffer);
  std::vector<RealType> runningratio(3,1.0);
  int sh(0);
  int nshells(ratios.size()/3);
  std::vector<RealType> orb_ratios;
//                 accumulate ratios on the way there
  for(int itz(0); itz<(iats.size()-1); itz++)
  {
    int iat=iats[itz];
    W.makeMove(iat,dR[itz]);
    RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
    runningratio[0]*=ratio;
    runningratio[1]*=orb_ratios[0];
    runningratio[2]*=ratio/orb_ratios[0];
    W.acceptMove(iat);
    Psi.acceptMove(W,iat);
    while(itz+1==rs[sh])
    {
      ratios[sh]=runningratio[0];
      ratios[nshells+sh]=runningratio[1];
      ratios[nshells*2+sh++]=runningratio[2];
    }
  }
  //we get to reject the last move
  {
    int iat=iats[iats.size()-1];
    W.makeMove(iat,dR[iats.size()-1]);
    RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
    runningratio[0]*=ratio;
    runningratio[1]*=orb_ratios[0];
    runningratio[2]*=ratio/orb_ratios[0];
    W.rejectMove(iat);
    Psi.rejectMove(iat);
    while(nshells*2+sh < ratios.size())
    {
      ratios[sh]=runningratio[0];
      ratios[nshells+sh]=runningratio[1];
      ratios[nshells*2+sh++]=runningratio[2];
    }
  }
//                 put the particles back
  for(int itz(0); itz<iats.size(); itz++)
    dR[itz]*=-1;
  for(int itz(iats.size()-2); itz>=0; itz--)
  {
    int iat=iats[itz];
    W.makeMove(iat,dR[itz]);
    RealType ratio = Psi.ratio(W,iat);
    W.acceptMove(iat);
    Psi.acceptMove(W,iat);
  }
//     Psi.updateBuffer(W,w_buffer,false);
//     W.saveWalker(thisWalker);
  return runningratio[0];
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

VMCUpdatePbyPWithDriftFast::RealType VMCUpdatePbyPWithDriftFast::advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios)
{
  Walker_t::Buffer_t& w_buffer(w1.DataSet);
  W.loadWalker(w1,true);
  Psi.copyFromBuffer(W,w_buffer);
  std::vector<RealType> logs;
  Psi.getLogs(logs);
  std::vector<RealType> runningratio(4,1.0);
  int sh(0);
//     runningratio[3]=std::exp(2.0*(logs[0]-csoffset));
  int nshells(rs.size());
  std::vector<RealType> orb_ratios;
//                 accumulate ratios on the way there
  for(int itz(0); itz<(iats.size()-1); itz++)
  {
    int iat=iats[itz];
    W.makeMove(iat,dR[itz]);
    RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
    runningratio[0]*=ratio;
    runningratio[1]*=orb_ratios[0];
    runningratio[2]*=orb_ratios[1];
    runningratio[3]*=ratio/orb_ratios[0];
    W.acceptMove(iat);
    Psi.acceptMove(W,iat);
    while(itz+1==rs[sh])
    {
      ratios[sh]=runningratio[0];
      ratios[nshells+sh]=runningratio[1];
      ratios[nshells*2+sh]=runningratio[2];
      ratios[nshells*3+sh++]=runningratio[3];
    }
  }
  //we get to reject the last move
  {
    int iat=iats[iats.size()-1];
    W.makeMove(iat,dR[iats.size()-1]);
    RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
    runningratio[0]*=ratio;
    runningratio[1]*=orb_ratios[0];
    runningratio[2]*=orb_ratios[1];
    runningratio[3]*=ratio/orb_ratios[0];
    W.rejectMove(iat);
    Psi.rejectMove(iat);
    while(nshells*3+sh < ratios.size())
    {
      ratios[sh]=runningratio[0];
      ratios[nshells+sh]=runningratio[1];
      ratios[nshells*2+sh]=runningratio[2];
      ratios[nshells*3+sh++]=runningratio[3];
    }
  }
//                 put the particles back
  for(int itz(0); itz<iats.size(); itz++)
    dR[itz]*=-1;
  for(int itz(iats.size()-2); itz>=0; itz--)
  {
    int iat=iats[itz];
    W.makeMove(iat,dR[itz]);
    RealType ratio = Psi.ratio(W,iat);
    W.acceptMove(iat);
    Psi.acceptMove(W,iat);
  }
//     Psi.updateBuffer(W,w_buffer,false);
//     W.saveWalker(thisWalker);
  return std::exp(2.0*(logs[0]-csoffset));
}


void VMCUpdatePbyPWithDriftFast::advanceWalker(Walker_t& thisWalker)
{
  Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
  W.loadWalker(thisWalker,true);
  Psi.copyFromBuffer(W,w_buffer);
  myTimers[1]->start();
  bool moved = false;
  for (int iter=0; iter<nSubSteps; ++iter)
  {
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    moved = false;
    for(int ig=0; ig<W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau*MassInvS[ig];
      RealType oneover2tau = 0.5/(tauovermass);
      RealType sqrttau = std::sqrt(tauovermass);
      for (int iat=W.first(ig); iat<W.last(ig); ++iat)
      {
        GradType grad_now=Psi.evalGrad(W,iat), grad_new;
        mPosType dr;
        getScaledDrift(tauovermass,grad_now,dr);
        dr += sqrttau*deltaR[iat];
        if (!W.makeMoveAndCheck(iat,dr))
        {
          ++nReject;
          continue;
        }
        //PosType newpos = W.makeMove(iat,dr);
        RealType ratio = Psi.ratioGrad(W,iat,grad_new);
        RealType prob = ratio*ratio;
        //zero is always rejected
        if (prob<std::numeric_limits<RealType>::epsilon())
        {
          ++nReject;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
          continue;
        }
        RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
        getScaledDrift(tauovermass,grad_new,dr);
        dr = thisWalker.R[iat]-W.R[iat]-dr;
        RealType logGb = -oneover2tau*dot(dr,dr);
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
    }
    //for subSteps must update thiswalker
    thisWalker.R=W.R;
    thisWalker.G=W.G;
    thisWalker.L=W.L;
  }
  myTimers[1]->stop();
  myTimers[2]->start();
  RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
  W.saveWalker(thisWalker);
  myTimers[2]->stop();
  myTimers[3]->start();
  EstimatorRealType eloc=H.evaluate(W);
  myTimers[3]->stop();
  //thisWalker.resetProperty(std::log(std::abs(psi)), psi,eloc);
  thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
  H.auxHevaluate(W,thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());
#if !defined(REMOVE_TRACEMANAGER)
  Traces->buffer_sample(W.current_step);
#endif
  if(!moved)
    ++nAllRejected;
}

void VMCUpdatePbyPWithDriftFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  myTimers[0]->start();
  for (; it != it_end; ++it)
  {
    advanceWalker(**it);
  }
  myTimers[0]->stop();
}


}

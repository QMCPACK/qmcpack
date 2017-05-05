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
    
    


#include <cmath>
#include <algorithm>
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Message/CommOperatorsMPI.h"
#include "OhmmsData/XMLcxx11Helper.h"
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


//   VMCUpdatePbyP::RealType VMCUpdatePbyP::advanceWalkerForCSEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios, std::vector<RealType>& weights, std::vector<RealType>& logs)
//   {
//     Walker_t::Buffer_t& w_buffer(w1.DataSet);
//     W.loadWalker(w1,true);
//     Psi.copyFromBuffer(W,w_buffer);
//     Psi.getLogs(logs);
//
//     RealType runningratio(std::exp(2.0*logs[1]));
//     int sh(0);
//     std::vector<RealType> lzratios;
// //                 accumulate ratios on the way there
//     for(int itz(0); itz<(iats.size()-1); itz++)
//     {
//       int iat=iats[itz];
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=lzratios[1];
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//       while(itz+1==rs[sh])
//         ratios[sh++]=runningratio;
//     }
//     //we get to reject the last move
//     {
//       int iat=iats[iats.size()-1];
//       W.makeMove(iat,dR[iats.size()-1]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=lzratios[1];
//       W.rejectMove(iat);
//       Psi.rejectMove(iat);
//       while(sh<ratios.size())
//         ratios[sh++]=runningratio;
//     }
// //                 put the particles back
//     for(int itz(0); itz<iats.size(); itz++)
//       dR[itz]*=-1;
//
//     for(int itz(iats.size()-2); itz>=0; itz--)
//     {
//       int iat=iats[itz];
//
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratio(W,iat);
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//     }
//
// //     Psi.updateBuffer(W,w_buffer,false);
// //     W.saveWalker(thisWalker);
//
//     return runningratio;
//   }

//   void VMCUpdatePbyP::advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone
//       , std::vector<MCWalkerConfiguration*>& wclone
//       , std::vector<QMCHamiltonian*>& hclone
//       , std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i)
//   {
//     int NumThreads(pclone.size());
//
//     //this can be modified for cache etc
//     RealType psi2_i_new[64];
//     for (int ipx=0; ipx<NumThreads; ++ipx)
//       psi2_i_new[ipx] = 2.0*W[ipx]->getPropertyBase()[LOGPSI] + c_i[ipx];
//
// #pragma omp parallel
//     {
//       RealType psi2_new=1.0;
//       for (int ipx=1; ipx<NumThreads; ++ipx)
//         psi2_new += std::exp(psi2_i_new[ipx]-psi2_i_new[0]);
//
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
//
//       //copy the new to now
//       RealType psi2_now=psi2_new;
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
//
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
//
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           mPosType dr=m_sqrttau*deltaR[iat];
//
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
// #pragma omp barrier
//           if (!movePtcl) continue;
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// // #pragma omp barrier
//           psi2_i_new[ip] = 2.0*std::log(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
//           psi2_new=1.0;
//           for(int jp=1; jp<NumThreads; ++jp)
//             psi2_new+=std::exp(psi2_i_new[jp]-psi2_i_new[0]);
//           RealType prob = std::exp(psi2_i_new[0]-psi2_0_now)*(psi2_new/psi2_now);
//           if (prob > rng_loc())
//           {
//             wclone[ip]->acceptMove(iat);
//             pclone[ip]->acceptMove(*wclone[ip],iat);
//
//             psi2_i_now=psi2_i_new[ip];
//             psi2_0_now=psi2_i_new[0];
//             psi2_now=psi2_new;
//           }
//           else
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//         }//loop over iat
//       }
//
//       Walker_t& thisWalker(*W[ip]);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       wclone[ip]->saveWalker(*W[ip]);
//       RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//       thisWalker.Weight=1.0;
//       hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
//
//     myTimers[0]->stop();
//   }

//   void VMCUpdatePbyP::estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//       , std::vector<MCWalkerConfiguration*>& wclone
//       , std::vector<QMCHamiltonian*>& hclone
//       , std::vector<RandomGenerator_t*>& rng
//       , std::vector<RealType>& ratio_i_0)
//   {
//     int NumThreads(pclone.size());
//
//     //this can be modified for cache etc
//     RealType psi2_i_new[64];
//     for (int ip=0; ip<NumThreads; ++ip)
//       psi2_i_new[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
//
//     RealType nn = -std::log(1.0*nSubSteps*W.getTotalNum());
// #pragma omp parallel
//     {
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
//
//       //copy the new to now
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
//
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
//
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           mPosType dr=m_sqrttau*deltaR[iat];
//
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
//           if (!movePtcl) continue;
//
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// #pragma omp barrier
//           psi2_i_new[ip] = 2.0*logl(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
// #pragma omp critical
//           {
//             ratio_i_0[ip] += expl( psi2_i_new[ip]-psi2_i_new[0] + nn);
//           }
//           wclone[ip]->rejectMove(iat);
//           pclone[ip]->rejectMove(iat);
//         }
//       }
//
//       //     Walker_t& thisWalker(*W[ip]);
//       //     Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       //     RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       //     wclone[ip]->saveWalker(*W[ip]);
//       //     RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //     //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       //     thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//
//       //     hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       //     hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
//
//     myTimers[0]->stop();
//   }

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
      for(int ig=0; ig<W.groups(); ++ig) //loop over species
      {
        RealType tauovermass = Tau*MassInvS[ig];
        RealType oneover2tau = 0.5/(tauovermass);
        RealType sqrttau = std::sqrt(tauovermass);
        for (int iat=W.first(ig); iat<W.last(ig); ++iat)
        {
          mPosType dr;
          getScaledDrift(tauovermass,W.G[iat],dr);
          dr += sqrttau*deltaR[iat];
          if (!W.makeMoveAndCheck(iat,dr))//reject illegal moves
          {
            ++nReject;
            continue;
          }
          RealType ratio = Psi.ratio(W,iat,dG,dL);
          RealType prob = ratio*ratio;
          //zero is always rejected
          if (prob<std::numeric_limits<RealType>::epsilon())
          {
            ++nReject;
            W.rejectMove(iat);
            Psi.rejectMove(iat);
            continue;
          }
          G = W.G+dG;
          RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
          getScaledDrift(tauovermass,G[iat],dr);
          dr = thisWalker.R[iat]-W.R[iat]-dr;
          RealType logGb = -oneover2tau*dot(dr,dr);
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
      //for subSteps must update thiswalker
      thisWalker.R=W.R;
      thisWalker.G=W.G;
      thisWalker.L=W.L;
    }
    myTimers[1]->stop();
    //if (moved)
    {
      myTimers[2]->start();
      RealType logpsi = Psi.evaluateLog(W,w_buffer);
      W.saveWalker(thisWalker);
      myTimers[2]->stop();
      myTimers[3]->start();
      EstimatorRealType eloc=H.evaluate(W);
      myTimers[3]->stop();
      //thisWalker.resetProperty(std::log(std::abs(psi)), psi,eloc);
      thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    if(!moved)
      ++nAllRejected;
    //else
    //{
    //  ++nAllRejected;
    //  H.rejectedMove(W,thisWalker);
    //}
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

//   VMCUpdatePbyPWithDriftFast::RealType VMCUpdatePbyPWithDriftFast::advanceWalkerForCSEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios, std::vector<RealType>& weights, std::vector<RealType>& logs )
//   {
//     Walker_t::Buffer_t& w_buffer(w1.DataSet);
//     W.loadWalker(w1,true);
//     Psi.copyFromBuffer(W,w_buffer);
//     Psi.getLogs(logs);
//     RealType logpsi = Psi.getLogPsi();
//     RealType logpsiCS = -(logpsi-logs[0]-csoffset);
//
//     RealType runningratio(1);
//     RealType runningweight(1);
// //     RealType runningweight(std::exp(2.0*logpsiCS));
//     logs[0]=runningweight;
//     int sh(0);
//     std::vector<RealType> lzratios;
// //                 accumulate ratios on the way there
//     for(int itz(0); itz<(iats.size()-1); itz++)
//     {
//       int iat=iats[itz];
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=ratio;
//       runningweight*=lzratios[0];
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//       while(itz+1==rs[sh])
//       {
//         ratios[sh]=runningratio;
//         weights[sh++]=runningweight;
//       }
//     }
//     //we get to reject the last move
//     {
//       int iat=iats[iats.size()-1];
//       W.makeMove(iat,dR[iats.size()-1]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=ratio;
//       runningweight*=lzratios[0];
//       W.rejectMove(iat);
//       Psi.rejectMove(iat);
//       while(sh<ratios.size())
//       {
//         ratios[sh]=runningratio;
//         weights[sh++]=runningweight;
//       }
//     }
// //                 put the particles back
//     for(int itz(0); itz<iats.size(); itz++)
//       dR[itz]*=-1;
//
//     for(int itz(iats.size()-2); itz>=0; itz--)
//     {
//       int iat=iats[itz];
//
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratio(W,iat);
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//     }
//
// //     Psi.updateBuffer(W,w_buffer,false);
// //     W.saveWalker(thisWalker);
//
//     return runningratio;
//   }

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
  myTimers[0]->stop();
}


/// Constructor.
VMCUpdateRenyiWithDriftFast::VMCUpdateRenyiWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg) :
  QMCUpdateBase(w,psi,h,rg)
{
  add_vmc_timers(myTimers);
}

VMCUpdateRenyiWithDriftFast::~VMCUpdateRenyiWithDriftFast()
{
}


void VMCUpdateRenyiWithDriftFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  myTimers[0]->start();
  WalkerIter_t begin(it);
  for (; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
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
    //Always compute the energy
    //if(moved)
    {
      myTimers[2]->start();
      //thisWalker.R = W.R;
      //w_buffer.rewind();
      //W.updateBuffer(w_buffer);
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
    }
    if(!moved)
      ++nAllRejected;
    //else
    //{
    //  ++nAllRejected;
    //  H.rejectedMove(W,thisWalker);
    //}
  }
  myTimers[0]->stop();
}


/// Constructor.
// VMCCSUpdatePbyPWithDriftFast::VMCCSUpdatePbyPWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//     QMCHamiltonian& h, RandomGenerator_t& rg) :
//     QMCUpdateBase(w,psi,h,rg)
// {
//   add_vmc_timers(myTimers);
// }
//
// VMCCSUpdatePbyPWithDriftFast::~VMCCSUpdatePbyPWithDriftFast()
// {
// }

//   void VMCUpdatePbyPWithDriftFast::advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i)
//   {
//     int NumThreads(pclone.size());
//     myTimers[0]->start();
//     myTimers[1]->start();
//     bool moved = false;
//
//     std::vector<RealType> psi2_i_now(NumThreads);
//     RealType psi2_now=0;
//     for (int ip=0; ip<NumThreads; ++ip)
//     {
//       psi2_i_now[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
//       psi2_now += std::exp(psi2_i_now[ip]-psi2_i_now[0]);
//     }
//
//     for (int iter=0; iter<nSubSteps; ++iter)
//     {
//       //create a 3N-Dimensional Gaussian with variance=1
//       makeGaussRandomWithEngine(deltaR, RandomGen);
//       moved = false;
//
//       for (int iat=0; iat<W.getTotalNum(); ++iat)
//       {
//         std::vector<GradType> grad_i_now(NumThreads);
//
//         // #pragma omp parallel for
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_i_now[ip] = pclone[ip]->evalGrad(*wclone[ip],iat);
//
//
//         GradType grad_now;
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_now += std::exp(psi2_i_now[ip]-psi2_i_now[0])/psi2_now*grad_i_now[ip];
//
//         mPosType dr;
//         getScaledDrift(m_tauovermass,grad_now,dr);
//         dr += m_sqrttau*deltaR[iat];
//
//
//         for (int ip=1; ip<NumThreads; ++ip)
//           wclone[ip]->makeMoveAndCheck(iat,dr);
//
//
//         if (!wclone[0]->makeMoveAndCheck(iat,dr))
//         {
//           ++nReject;
//           continue;
//         }
//
//         RealType psi2_new(0);
//         std::vector<GradType> grad_i_new(NumThreads);
//         std::vector<RealType> psi2_i_new(NumThreads);
//         // #pragma omp parallel for
//         for (int ip=0; ip<NumThreads; ++ip)
//         {
//           RealType ratio = pclone[ip]->ratioGrad(*wclone[ip],iat,grad_i_new[ip]);
//           psi2_i_new[ip] = 2.0*std::log(std::abs(ratio)) + psi2_i_now[ip];
//         }
//         for (int ip=0; ip<NumThreads; ++ip) psi2_new += std::exp(psi2_i_new[ip]-psi2_i_new[0]);
//
//         RealType prob = std::exp(psi2_i_new[0]-psi2_i_now[0])*(psi2_new/psi2_now);
//
//         //zero is always rejected
//         if (prob < std::numeric_limits<RealType>::epsilon())
//         {
//           ++nReject;
//           for (int ip=0; ip<NumThreads; ++ip)
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//           continue;
//         }
//
//         RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
//
//         GradType grad_new;
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_new += std::exp(psi2_i_new[ip]-psi2_i_new[0])/psi2_new*grad_i_new[ip];
//
//         getScaledDrift(m_tauovermass,grad_new,dr);
//         dr += W[0]->R[iat] - wclone[0]->R[iat] -dr;
//
//         RealType logGb = -m_oneover2tau*dot(dr,dr);
//
//         if ( RandomGen() < prob*std::exp(logGb-logGf))
//         {
//           moved = true;
//           ++nAccept;
//           for (int ip=0; ip<NumThreads; ++ip)
//           {
//             wclone[ip]->acceptMove(iat);
//             pclone[ip]->acceptMove(*wclone[ip],iat);
//             psi2_i_now[ip]=psi2_i_new[ip];
//           }
//           psi2_now = psi2_new;
//         }
//         else
//         {
//           ++nReject;
//           for (int ip=0; ip<NumThreads; ++ip)
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//         }
//       }
//       //for subSteps must update walkers
//       for (int ip=0; ip<NumThreads; ++ip) wclone[ip]->saveWalker(*W[ip]);
//     }
//
//     myTimers[1]->stop();
//
//     // #pragma omp parallel for
//     for (int ip=0; ip<NumThreads; ++ip)
//     {
//       Walker_t& thisWalker(*W[ip]);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       if (moved)
//       {
//         RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//         wclone[ip]->saveWalker(*W[ip]);
//         EstimatorRealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//         //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//         thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//         hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//         hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//       }
//       else
//       {
//         ++nAllRejected;
//         hclone[ip]->rejectedMove(*wclone[ip],thisWalker);
//       }
//     }
//
//     myTimers[0]->stop();
//   }

//   VMCUpdatePbyPSampleRN::VMCUpdatePbyPSampleRN(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide,
//       QMCHamiltonian& h, RandomGenerator_t& rg):
//     QMCUpdateBase(w,psi,guide,h,rg), logEpsilon(0.0)
//   {
//     add_vmc_timers(myTimers);
//     //   dummy= new MCWalkerConfiguration(w);
//   }
//
//   VMCUpdatePbyPSampleRN::~VMCUpdatePbyPSampleRN()
//   {
//   }
//
//   void VMCUpdatePbyPSampleRN::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
//   {
//     UpdatePbyP=true;
//     //   Guide.resizeTempP(W);
//
//     for (;it != it_end; ++it)
//     {
//       Walker_t& thisWalker(**it);
//       W.loadWalker(thisWalker,UpdatePbyP);
//
//       Walker_t::Buffer_t tbuffer;
//       RealType logguide=Guide.registerData(W,tbuffer)+logEpsilon;
//       RealType logpsi=Psi.registerData(W,tbuffer);
//       thisWalker.DataSet=tbuffer;
//       RealType ene = H.evaluate(W);
//
//       thisWalker.resetProperty(logpsi,Psi.getPhase(),ene);
//       H.saveProperty(thisWalker.getPropertyBase());
//       thisWalker.ReleasedNodeAge=0;
//       thisWalker.ReleasedNodeWeight=0;
//       //       thisWalker.Weight=1.0;
//       thisWalker.Weight=std::exp(2.0*(logpsi-logguide));
//     }
//   }
//
//   void VMCUpdatePbyPSampleRN::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
//
//     myTimers[0]->start();
//     for (; it != it_end; ++it)
//     {
//
//       Walker_t& thisWalker(**it);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//
//       W.loadWalker(thisWalker,true);
//       //       dummy.loadWalker(thisWalker,true);
//       //       RealType logguide2_now=2.0*(logEpsilon+Guide.evaluateLog(W));
//       Guide.copyFromBuffer(W,w_buffer);
//       RealType logguide2_now = 2.0*(logEpsilon+Guide.getLogPsi());
//
//       Psi.copyFromBuffer(W,w_buffer);
//       RealType logpsi2_now = 2.0*Psi.getLogPsi();
//
//       myTimers[1]->start();
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         makeGaussRandomWithEngine(deltaR,RandomGen);
//         bool stucked=true;
//         for (int iat=0; iat<W.getTotalNum(); ++iat)
//         {
//           mPosType dr = m_sqrttau*deltaR[iat];
//           //ignore illegal moves
//           if (!W.makeMoveAndCheck(iat,dr))
//           {
//             ++nReject;
//             continue;
//           }
//           //               dummy.makeMoveAndCheck(iat,dr);
//           //PosType newpos = W.makeMove(iat,dr);
//           RealType psi_ratio = Psi.ratio(W,iat);
//           RealType logpsi2_new = 2.0*std::log(std::abs(psi_ratio))+logpsi2_now;
//           RealType logguide2_new = 2.0*std::log(Guide.ratio(W,iat))+ logguide2_now;
//
//           long double prob = psi_ratio*psi_ratio*(1.0 + expl(logguide2_new-logpsi2_new))/(1.0+ expl(logguide2_now-logpsi2_now));
//           //               app_log()<<prob<< std::endl;
//           //RealType prob = std::min(1.0e0,ratio*ratio);
//           if (RandomGen() < prob)
//           {
//             stucked=false;
//             ++nAccept;
//             W.acceptMove(iat);
//             Guide.acceptMove(W,iat);
//             Psi.acceptMove(W,iat);
//             logpsi2_now=logpsi2_new;
//             logguide2_now=logguide2_new;
//           }
//           else
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//             Guide.rejectMove(iat);
//           }
//         }
//         if (stucked)
//         {
//           ++nAllRejected;
//           H.rejectedMove(W,thisWalker);
//         }
//       }
//       myTimers[1]->stop();
//
//       myTimers[2]->start();
//       //thisWalker.R = W.R;
//       //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
//       //w_buffer.rewind();
//       //W.updateBuffer(w_buffer);
//       //       RealType logguide = logguide2_now;
//       RealType logguide = Guide.updateBuffer(W,w_buffer,true);
//       RealType logpsi = Psi.updateBuffer(W,w_buffer,true);
//       W.saveWalker(thisWalker);
//
//       myTimers[2]->stop();
//
//       myTimers[3]->start();
//       EstimatorRealType eloc=H.evaluate(W);
//       myTimers[3]->stop();
//
//       thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
//       //       thisWalker.resetProperty(0.5*logpsi2_now,Psi.getPhase(),eloc);
//       H.auxHevaluate(W,thisWalker);
//       //       thisWalker.Weight = 1.0;
//       //       thisWalker.Weight = expl(logpsi2_now-logguide2_now);
//       thisWalker.Weight = 1.0/(1.0+expl(2.0*(logguide+logEpsilon-logpsi)));
//       H.saveProperty(thisWalker.getPropertyBase());
//     }
//     myTimers[0]->stop();
//   }
//
//   void VMCUpdatePbyPSampleRN::advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone
//       , std::vector<MCWalkerConfiguration*>& wclone
//       , std::vector<QMCHamiltonian*>& hclone
//       , std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i)
//   {
//     int NumThreads(pclone.size());
//
//     //this can be modified for cache etc
//     long double psi2_i_new[64];
//     for (int ip=0; ip<NumThreads; ++ip)
//       psi2_i_new[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI] + c_i[ip];
//
// #pragma omp parallel
//     {
//       long double psi2_new=1.0;
//       for (int ip=1; ip<NumThreads; ++ip)
//         psi2_new += std::exp(psi2_i_new[ip]-psi2_i_new[0]);
//
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
//
//       //copy the new to now
//       long double psi2_now=psi2_new;
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
//
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
//
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           mPosType dr=m_sqrttau*deltaR[iat];
//
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
//           if (!movePtcl) continue;
//
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// #pragma omp barrier
//           psi2_i_new[ip] = 2.0*std::log(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
//           psi2_new=1.0;
//           for(int jp=1; jp<NumThreads; ++jp)
//             psi2_new+= expl(psi2_i_new[jp]-psi2_i_new[0]);
//           long double prob = expl(psi2_i_new[0]-psi2_0_now)*(psi2_new/psi2_now);
//           if (prob > rng_loc())
//           {
//             wclone[ip]->acceptMove(iat);
//             pclone[ip]->acceptMove(*wclone[ip],iat);
//
//             psi2_i_now=psi2_i_new[ip];
//             psi2_0_now=psi2_i_new[0];
//             psi2_now=psi2_new;
//           }
//           else
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//         }//loop over iat
//       }
//
//       Walker_t& thisWalker(*W[ip]);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       wclone[ip]->saveWalker(*W[ip]);
//       RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//       thisWalker.Weight=1.0;
//       hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
//
//     myTimers[0]->stop();
//   }
//
//   void VMCUpdatePbyPSampleRN::estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//       , std::vector<MCWalkerConfiguration*>& wclone
//       , std::vector<QMCHamiltonian*>& hclone
//       , std::vector<RandomGenerator_t*>& rng
//       , std::vector<RealType>& ratio_i_0)
//   {
//     int NumThreads(pclone.size());
//
//     //this can be modified for cache etc
//     RealType psi2_i_new[64];
//     for (int ip=0; ip<NumThreads; ++ip)
//       psi2_i_new[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
//
// //     RealType nn = -std::log(1.0*nSubSteps*W.getTotalNum());
// #pragma omp parallel
//     {
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
//
//       //copy the new to now
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
//
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
//
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           mPosType dr=m_sqrttau*deltaR[iat];
//
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
//           if (!movePtcl) continue;
//
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// #pragma omp barrier
//           psi2_i_new[ip] = 2.0*logl(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
// #pragma omp critical
//           {
//             ratio_i_0[ip] += expl( psi2_i_new[ip]-psi2_i_new[0]);
//           }
//           wclone[ip]->rejectMove(iat);
//           pclone[ip]->rejectMove(iat);
//         }
//       }
//
//       //     Walker_t& thisWalker(*W[ip]);
//       //     Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       //     RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       //     wclone[ip]->saveWalker(*W[ip]);
//       //     RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //     //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       //     thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//
//       //     hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       //     hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
//     RealType nn = 1.0/ratio_i_0[0];
//     for (int ip=0; ip<NumThreads; ++ip) ratio_i_0[ip]*=nn;
//     myTimers[0]->stop();
//   }

// initialize static data members
std::vector<std::vector<VMCUpdatePbyPNodeless::RealType> > VMCUpdatePbyPNodeless::tfl_history;
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::tfl_avg = 0.0;
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::tfl_sdv = -1.0;
//bool VMCUpdatePbyPNodeless::setNGMag = false;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Constructor for the VMCUpdate class that uses particle-by-particle moves and a
///         nodeless guiding function
///
/// \param[in,out]  w      the set of walkers to use
/// \param[in,out]  psi    the trial function, which in general has nodes
/// \param[in,out]  h      the Hamiltonian
/// \param[in,out]  rg     random number generator
/// \param[in]      ips    ion particle set (tells us where the ions are positioned)
/// \param[in]      eps    coefficient for nodeless guiding adjustment
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless(MCWalkerConfiguration & w,
                                             TrialWaveFunction & psi,
                                             QMCHamiltonian & h,
                                             RandomGenerator_t & rg,
                                             const ParticleSet & ips,
                                             const RealType eps)
: QMCUpdateBase(w,psi,h,rg)
, NodelessEpsilon(eps)
, IonPositions(ips.R.begin(), ips.R.end())
{

  // ensure some sanity
  if ( NodelessEpsilon < RealType(0) )
    APP_ABORT("VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless was given a negative epsilon");
  //if ( NodelessAlpha < RealType(0) )
  //  APP_ABORT("VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless was given a negative alpha");
  //if ( NodelessBeta < RealType(0) )
  //  APP_ABORT("VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless was given a negative beta");

  //// print the ion positions
  //std::cout << std::endl;
  //std::cout << "printing ion positions in VMCUpdatePbyPNodeless::VMCUpdatePbyPNodeless" << std::endl;
  //for (auto it = IonPositions.begin(); it != IonPositions.end(); it++)
  //  std::cout << *it << std::endl;
  //std::cout << std::endl;
  //APP_ABORT("Stopping after printing ion particle positions for Nodeless guiding");

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Read parameters of the nodeless guiding function from xml
///
/// \param[in,out]  cur      pointer to the xml node to read from
///
///////////////////////////////////////////////////////////////////////////////////////////////////
bool VMCUpdatePbyPNodeless::put(xmlNodePtr cur) {

  // do the usual base class put
  const bool base_success = QMCUpdateBase::put(cur);

  // find and process nodeless guiding node
  bool my_success = false;
  for ( auto kur : getXMLChildList(cur, "nodelessGuiding") ) {

    // remember that we found something
    my_success = true;

    // reset guiding function parameters
    cgCountSigmas.clear();
    cgCountNelecs.clear();
    cgGaussStarts.clear();
    cgGaussEnds.clear();
    cgGaussAlphas.clear();
    cgGaussSigmas.clear();
    cgGaussCenters.clear();
    mdBetas.clear();
    mdDists.clear();
    mdCenters.clear();

    // read counting groups' info
    for ( auto chl : getXMLChildList(kur, "countGroup") ) {

      // read counting group standard deviation and target number of electrons
      cgCountSigmas.push_back(-1.0);
      cgCountNelecs.push_back(0.0);
      getXMLAttributes(chl, *cgCountSigmas.rbegin(), "sigma", *cgCountNelecs.rbegin(), "nelec");

      // record the start of this counting group's gaussians
      cgGaussStarts.push_back(cgGaussSigmas.size());

      // read the standard deviations and centers of the gaussians
      for ( auto khl : getXMLChildList(chl, "gaussian") ) {

        // <gaussian alpha="1.0" sigma="1.3" type="Array"> 0.1 0.2 0.3 </gaussian>

        // read alpha and sigma values for this gaussian
        cgGaussAlphas.push_back(-1.0);
        cgGaussSigmas.push_back(-1.0);
        getXMLAttributes(khl, *cgGaussAlphas.rbegin(), "alpha", *cgGaussSigmas.rbegin(), "sigma");

        // read the coordinates of the gaussian's center
        cgGaussCenters.push_back(ParticleSet::SingleParticlePos_t());
        if ( ! putContent(cgGaussCenters.rbegin()->begin(), cgGaussCenters.rbegin()->end(), khl) )
          APP_ABORT("ERROR: problem reading coordinates of a countGroup gaussian center");

      }

      // record the end of this counting group's gaussians
      cgGaussEnds.push_back(cgGaussSigmas.size());

    }

    // read minimum distance information
    for ( auto chl : getXMLChildList(kur, "minDistCenter") ) {

      // <minDistCenter beta="1.0" dist="1.3"> -3.1 3.3 4.7 </minDistCenter>

      // read switching speed and distance parameters
      mdBetas.push_back(-1.0);
      mdDists.push_back(-1.0);
      getXMLAttributes(chl, *mdBetas.rbegin(), "beta", *mdDists.rbegin(), "dist");

      // read the coordinates of the center to measure distance from
      mdCenters.push_back(ParticleSet::SingleParticlePos_t());
      if ( ! putContent(mdCenters.rbegin()->begin(), mdCenters.rbegin()->end(), chl) )
        APP_ABORT("ERROR: problem reading coordinates of a min distance center");

    }

    // check that things that should be positive are
    check_all_positive(cgCountSigmas, "cgCountSigmas");
    check_all_positive(cgGaussAlphas, "cgGaussAlphas");
    check_all_positive(cgGaussSigmas, "cgGaussSigmas");
    check_all_positive(mdBetas, "mdBetas");
    check_all_positive(mdDists, "mdDists");

    // check that there is at least one minimum distance center
    if ( mdBetas.empty() )
      APP_ABORT("ERROR: Found no minimum distance centers while reading nodelessGuiding entry in xml input");

    // print what we have
    if ( OHMMS::Controller->rank() == 0 && omp_get_thread_num() == 0 ) {
      app_log() << std::endl;
      app_log() << omp_get_num_threads() << " threads are running through VMCUpdatePbyPNodeless::put" << std::endl;
      app_log() << std::endl;
      app_log() << "NodelessEpsilon = " << NodelessEpsilon << std::endl;
      app_log() << std::endl;
      for (int i = 0; i < cgCountSigmas.size(); i++) {
        app_log() << "Count group:" << std::endl;
        app_log() << "  count sigma = " << cgCountSigmas.at(i) << std::endl;
        app_log() << "  count nelec = " << cgCountNelecs.at(i) << std::endl;
        app_log() << "  Gaussians" << std::endl;
        for (int j = cgGaussStarts.at(i); j < cgGaussEnds.at(i); j++) {
          app_log() << "    alpha = " << cgGaussAlphas.at(j) << std::endl;
          app_log() << "    sigma = " << cgGaussSigmas.at(j) << std::endl;
          app_log() << "   center = " << cgGaussCenters.at(j) << std::endl;
        }
        app_log() << std::endl;
      }
      for (int i = 0; i < mdBetas.size(); i++) {
        app_log() << "Min Dist Center:" << std::endl;
        app_log() << "    beta = " << mdBetas.at(i) << std::endl;
        app_log() << "    dist = " << mdDists.at(i) << std::endl;
        app_log() << "  center = " << mdCenters.at(i) << std::endl;
        app_log() << std::endl;
      }
    }

  }

  // make sure we found something
  if (!my_success)
    APP_ABORT("ERROR: failed to find nodelessGuiding entry in xml input");

  // return whether success was had
  return base_success && my_success;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Empties the history vectors and ensures they are the correct size.
///
/// \param[in,out]  comm     communicator for averaging over all processes
/// \param[in]      nthread  number of threads per process
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::reset_history(Communicate * const comm, const int nthread)
{

  // ensure there is one history vector per thread
  if ( tfl_history.size() < nthread )
    tfl_history.resize(nthread);

  // empty the history vectors
  for (int i = 0; i < tfl_history.size(); i++)
    tfl_history.at(i).clear();

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes and stores the average and standard deviations of the history of trial
///         funciton logarithms and then clears the history.
///
/// \param[in,out]  comm     communicator for averaging over all processes
/// \param[in]      nthread  number of threads per process
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::process_history(Communicate * const comm, const int nthread)
{

  // compute the average and std dev of the trial function logpsi
  // values that will be used in nodeless guiding
  std::vector<RealType> sums(3, 0.0);
  for (int i = 0; i < tfl_history.size(); i++) {
    for (int j = 0; j < tfl_history[i].size(); j++) {
      const RealType tfl = tfl_history[i][j];
      sums[0] += tfl;
      sums[1] += tfl * tfl;
      sums[2] += 1.0;
    }
  }
  comm->allreduce(sums);
  if ( sums[2] > 0.1 ) // do nothing if we have no previous sample data to work with
  {

    const RealType temp_avg = sums[0] / sums[2];
    const RealType temp_sdv = std::sqrt( sums[1] / sums[2] - temp_avg * temp_avg );
    app_log() << std::endl;
    app_log() << "Estimating average of trial function logarithm for nodeless guiding: " << std::endl;
    app_log() << "  samples = " << sums[2] << std::endl;
    app_log() << "  tfl_avg = " << temp_avg << std::endl;
    app_log() << "  tfl_sdv = " << temp_sdv << std::endl;
    app_log() << std::endl;

    // only record the values if we don't have values already
    if ( tfl_sdv < 0.0 ) {
      tfl_avg = temp_avg;
      tfl_sdv = temp_sdv;
      app_log() << "Nodeless guiding magnitude set based on this tfl_avg" << std::endl;
      app_log() << std::endl;
    }

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Resets info stored for trial function logarithm.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::reset_tfl()
{
  tfl_avg = 0.0;
  tfl_sdv = -1.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Destructor for the VMCUpdate class that uses particle-by-particle moves and a
///         nodeless guiding function
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::~VMCUpdatePbyPNodeless() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Initialize internal nodeless guiding info for the provided configuration.
///
/// \param[in]      P        holds the configuration information
/// \param[in]      tfl      trial function logarithm
///
/// \return  the value of the overall guiding function after initialization
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::init_nodeless(const ParticleSet & P, const RealType tfl) {

  // get dimensions
  const int np = P.R.size(); // number of particles
  const int nc = cgCountSigmas.size(); // number of counting groups
  const int nd = mdBetas.size(); // number of min distance centers

  // get a temporary particle position that will be useful
  ParticleSet::SingleParticlePos_t tpp;

  // initialize un-normalized counting functions and their normalizing constants
  cgUnormalized.assign(np * nc, 0);
  cgNorms.assign(np, 0);
  for (int i = 0; i < np; i++) {
    for (int k = 0; k < nc; k++) {
      for (int l = cgGaussStarts.at(k); l < cgGaussEnds.at(k); l++) {
        tpp = cgGaussCenters.at(l) - P.R[i];
        cgUnormalized.at(k+i*nc) += cgGaussAlphas.at(l) * std::exp( -dot(tpp,tpp) / ( 2.0 * cgGaussSigmas.at(l) * cgGaussSigmas.at(l) ) );
      }
      cgNorms.at(i) += cgUnormalized.at(k+i*nc);
    }
  }

  // compute each counting group's total count
  cgCounts.assign(nc, 0);
  for (int i = 0; i < np; i++)
    for (int k = 0; k < nc; k++)
      cgCounts.at(k) += cgUnormalized.at(k+i*nc) / cgNorms.at(i);

  // get sum of counting group penalty exponents
  cgPenaltyExponent = 0;
  for (int k = 0; k < nc; k++)
    cgPenaltyExponent -=   ( cgCountNelecs.at(k) - cgCounts.at(k) ) * ( cgCountNelecs.at(k) - cgCounts.at(k) )
                         / ( 2.0 * cgCountSigmas.at(k) * cgCountSigmas.at(k) );

  // get product of min distance penalties
  mdPenalty = 1;
  for (int i = 0; i < np; i++) {
    RealType max_val = 0;
    for (int k = 0; k < nd; k++) {
      tpp = mdCenters.at(k) - P.R[i];
      max_val = std::max(max_val, 1.0 / ( 1.0 + std::exp( mdBetas.at(k) * ( std::sqrt( std::abs( dot(tpp,tpp) ) ) - mdDists.at(k) ) ) ) );
    }
    mdPenalty *= max_val;
  }

  // initialize the nodeless adjustment as the epsilon-scaled "average" trial function value
  RealType nodelessAdj = NodelessEpsilon * std::exp( 2.0 * tfl_avg );

  // penalize the nodeless adjustment by the min distance and counting group penalties
  nodelessAdj *= mdPenalty * std::exp(cgPenaltyExponent);

  // If we don't have an average and standard deviation for a previously-taken set of
  // trial function logarithms, return the un-adjusted square norm of the trial function
  if ( tfl_sdv < 0.0 )
    return std::exp( 2.0 * tfl );

  // otherwise, return the trial function square norm plus the penalized nodeless adjustment
  return std::exp( 2.0 * tfl ) + nodelessAdj;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Update internal nodeless guiding info after a single particle move
///
/// \param[in]      P        holds the configuration information
/// \param[in]      iat      index of the moved particle
/// \param[in]      tfl      trial function logarithm after the move
///
/// \return  the value of the overall guiding function after the update
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::update_nodeless(const ParticleSet & P, const int iat, const RealType tfl) {
  return this->init_nodeless(P, tfl);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that moves the supplied group of walkers according to a Metropolis-Hastings
///         random walk over the nodeless guiding function.
///
/// \param[in,out]  it       iterator to the first walker to be moved
/// \param[in,out]  it_end   iterator to one past the last walker to be moved
/// \param[in]      measure  ??? (not sure what this is for, it does not get used here)
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void VMCUpdatePbyPNodeless::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{

  // initialize a counter for the number of walkers
  int nw = 0;

  // loop over the walkers
  for (; it != it_end; it++) {

    // get a reference to the current walker
    Walker_t& thisWalker(**it);

    // load the walker information into our particle set
    W.loadWalker(thisWalker,true);

    // get a reference to where this walker stores its trial function state
    Walker_t::Buffer_t & w_buffer(thisWalker.DataSet);

    // load trial wave function state for this walker into our trial wave function object
    Psi.copyFromBuffer(W, w_buffer);

    // initialize the old configuration's trial funciton logpsi and nodeless guiding function square norm
    Psi.evaluateLog(W);
    RealType old_tfl = Psi.getLogPsi();
    //RealType old_sqn = this->get_nodeless_gf_sqn(old_tfl, W.G, W.L);
    //RealType old_sqn = this->get_nodeless_gf_sqn_new(old_tfl);
    RealType old_sqn = this->init_nodeless(W, old_tfl);
    //std::cout << std::endl;

    // loop over the sub-steps in the sampling walk
    for (int si = 0; si < nSubSteps; si++) {

      // create this sub-step's random displacements for each particle 
      makeGaussRandomWithEngine(deltaR, RandomGen);

      // initialize flag to tell if any particles moved during this substep
      bool stuck = true;

      // loop over types of particles
      for(int ig = 0; ig < W.groups(); ig++) {

        // get the mass-modified time step
        RealType sqrttau = std::sqrt(Tau*MassInvS[ig]);

        // loop over particles of this species
        for (int iat = W.first(ig); iat < W.last(ig); iat++) {

          // get the proposed displacement for this particle
          const mPosType dr = sqrttau*deltaR[iat];

          // set up the move in our particle set
          const bool move_is_legal = W.makeMoveAndCheck(iat, dr);

          // reject illegal moves
          if ( ! move_is_legal ) {
            nReject++;
            continue;
          }

//          // get the trial function ratio and the change in the gradient and laplacian for the proposed move
//          const RealType ratio = Psi.ratio(W, iat, dG, dL);
//
//          // get the gradient and laplacian at the proposed configuration
//          G = W.G + dG;
//          L = W.L + dL;

          // get the trial function ratio
          const RealType ratio = Psi.ratio(W,iat);

          // get the log of the trial function at the new configuration
          const RealType new_tfl = std::log(std::abs(ratio)) + old_tfl;

          // get the square norm of the nodeless guiding function
          //const RealType new_sqn = this->get_nodeless_gf_sqn(new_tfl, G, L);
          //const RealType new_sqn = this->get_nodeless_gf_sqn_new(new_tfl);

          // update internal parameters to reflect the move and get the new square norm of the guiding function 
          const RealType new_sqn = this->update_nodeless(W, iat, new_tfl);

          // if the ratio of square norms satisfies the Metropolis condition, accept the move
          if ( RandomGen() < new_sqn / old_sqn ) {

            nAccept++;
            stuck = false;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
            //W.G = G;
            //W.L = L;
            old_tfl = new_tfl;
            old_sqn = new_sqn;

//            {
//              // get the sum of the laplacian and square gradient of log(Psi)
//              RealType d = 0.0;
//              for (int i = 0; i < W.G.size(); i++) {
//                d += dot(W.G[i], W.G[i]);
//                d += W.L[i];
//              }
//
//              // get |Psi|^2 and |laplacian of Psi|^2
//              const RealType psi2 = std::exp( 2.0 * new_tfl );
//              const RealType lap2 = psi2 * d * d;
//
//              std::printf("    PsiDel2Psi/g = %10.2e", std::sqrt(lap2*psi2) / new_sqn);
//            }

          // otherwise, reject the move
          } else {

            nReject++;
            W.rejectMove(iat);
            Psi.rejectMove(iat);

            // undo the update
            this->update_nodeless(W, iat, old_tfl);

          }

          //std::cout << std::endl;

        } // end iat loop over particles

      } // end ig loop over types of particles

      // if nothing moved, increment the all-rejected counter
      nAllRejected += ( stuck ? 1 : 0 );

      //// record our new position, gradient, and laplacian
      //thisWalker.R = W.R;
      //thisWalker.G = W.G;
      //thisWalker.L = W.L;

      // record our new position
      thisWalker.R = W.R;

    } // end si loop over sub-steps

    // for now, we are assuming there is only one walker to advance
    if ( ++nw > 1 )
      APP_ABORT("VMCUpdatePbyPNodeless::advanceWalkers encountered more than one walker");

    // save a history of the logarithm of the underlying trial function
    tfl_history.at(omp_get_thread_num()).push_back(old_tfl);

    // save the trial wave function state
    const RealType logpsi = Psi.updateBuffer(W, w_buffer, false);

    // save the particle set information
    W.saveWalker(thisWalker);

    // evaluate the local energy
    EstimatorRealType eloc = H.evaluate(W);

    //std::printf("psi^2 = %10.2e    g = %10.2e    eloc = %10.2e    toAvg = %10.2e\n", std::exp(2.0*logpsi), old_sqn, eloc, std::exp(2.0*logpsi) * eloc / old_sqn);

    // Save some basic info about the underlying trial function.
    // Note that here we save the logpsi value of the nodeless guiding function, NOT of the underlying trial function.
    thisWalker.resetProperty(0.5*std::log(old_sqn), Psi.getPhase(), eloc);

    // not sure what these do
    H.auxHevaluate(W, thisWalker);
    H.saveProperty(thisWalker.getPropertyBase());

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not implemented...
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::advanceWalkerForEE(Walker_t& w1,
                                                                          std::vector<PosType>& dR,
                                                                          std::vector<int>& iats,
                                                                          std::vector<int>& rs,
                                                                          std::vector<RealType>& ratios)
{
  APP_ABORT("VMCUpdatePbyPNodeless::advanceWalkerForEE not implemented");
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns the square norm of the nodeless guiding function.
///         This guiding function is defined to be
///
///         g = sqrt( psi^2 + epsilon * (laplacian of psi)^2 / ( 1 + exp( alpha * ( logpsi - mu + sigma ) / sigma ) ) )
///
///         where epsilon and alpha are user-chosen parameters and mu and sigma are the average and standard deviation
///         of a previously-sampled set of trial function logarithm values.
///
/// \param[in]      logpsi   the log of the trial function
/// \param[in]      grad     the gradient of the log of the trial function
/// \param[in]      lap      the laplacian of the log of the trial function
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::get_nodeless_gf_sqn(const RealType logpsi,
                                                                           const ParticleSet::ParticleGradient_t & grad,
                                                                           const ParticleSet::ParticleLaplacian_t & lap) const
{

  // get the sum of the laplacian and square gradient of log(Psi)
  RealType d = 0.0;
  for (int i = 0; i < grad.size(); i++) {
    d += dot(grad[i], grad[i]);
    d += lap[i];
  }

  // get |Psi|^2 and |laplacian of Psi|^2
  const RealType psi2 = std::exp( 2.0 * logpsi );
  const RealType lap2 = psi2 * d * d;

  // If we don't have an average and standard deviation for a previously-taken set of
  // trial function logarithms, use a simpler nodeless guiding function instead.
  // The idea is to use this for the very first warmup when we have no history
  // to work with and then to switch to the general formula after that.
  if ( tfl_sdv < 0.0 )
    return psi2 + NodelessEpsilon * lap2 / 100.0 ;

  // return the square norm of the nodeless guiding function
  return psi2 + NodelessEpsilon * lap2 / ( 1.0 + std::exp( NodelessAlpha * ( logpsi - tfl_avg + tfl_sdv ) / tfl_sdv ) );

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns the new nodeless guiding function.
///         This guiding function is defined to be
///
///                   psi^2 + epsilon * <psi^2> * prod_i ( 1 / ( 1 + exp( beta * ( rm_i - alpha ) ) ) )
///
///         where epsilon, beta, and alpha are user-chosen parameters rm_i is the distance from
///         the ith electron to the nearest nucleus
///
/// \param[in]      logpsi   the log of the trial function
///
///////////////////////////////////////////////////////////////////////////////////////////////////
VMCUpdatePbyPNodeless::RealType VMCUpdatePbyPNodeless::get_nodeless_gf_sqn_new(const RealType logpsi) const
{

  // If we don't have an average and standard deviation for a previously-taken set of
  // trial function logarithms, use a simpler nodeless guiding function instead.
  // The idea is to use this for the very first warmup when we have no history
  // to work with and then to switch to the general formula after that.
  if ( tfl_sdv < 0.0 )
    return std::exp( 2.0 * logpsi );

  // !!! THIS IS NOT QUITE RIGHT.  YOU NEED TO ACTUALLY AVERAGE PSI^2 !!!
  // initialize return based on average psi2
  RealType retval = NodelessEpsilon * std::exp( 2.0 * tfl_avg );

  ParticleSet::SingleParticlePos_t temp_vec;

  // loop over types of particles
  for(int ig = 0; ig < W.groups(); ig++) {

    //std::cout << "ig = " << ig << std::endl;

    // loop over particles of this species
    for (int iat = W.first(ig); iat < W.last(ig); iat++) {

      //auto temp_vec = 1.0 * W.R[iat];

      //std::cout << "  iat = " << iat;
      //std::cout << "    W.R[iat] = " << W.R[iat] << std::endl;

      // get the smallest electron-ion distance for this electron
      RealType min_dist = 1.0e100;
      for (auto it = IonPositions.begin(); it != IonPositions.end(); it++) {
        temp_vec = *it - W.R[iat];
        min_dist = std::min(min_dist, std::abs(std::sqrt(dot(temp_vec, temp_vec))));
      }

      //std::printf("  %10.2e", min_dist);

      // apply the distance penalty for this electron
      retval *= 1.0 / ( 1.0 + std::exp( NodelessBeta * ( min_dist - NodelessAlpha ) ) );

    }

  }

  // print what happened
  //std::printf("    psi^2 = %10.2e    adjust by %10.2e", std::exp( 2.0 * logpsi ), retval);
  //std::cout << std::endl;

  // add psi2
  retval += std::exp( 2.0 * logpsi );

  //APP_ABORT("VMCUpdatePbyPNodeless::get_nodeless_gf_sqn_new stopping here");

  // return nodeless guiding function value
  return retval;

}

}


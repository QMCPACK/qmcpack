//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateReny.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

namespace qmcplusplus
{

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

}

/***************************************************************************
 * $RCSfile: VMCUpdatePbyP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdatePbyP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/

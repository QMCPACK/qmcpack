//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

VMCUpdateAll::VMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                           QMCHamiltonian& h, RandomGenerator_t& rg)
  : QMCUpdateBase(w,psi,h,rg)
{
  UpdatePbyP=false;
}

VMCUpdateAll::~VMCUpdateAll()
{
}

void VMCUpdateAll::advanceWalker(Walker_t& thisWalker, bool recompute)
{
    makeGaussRandomWithEngine(deltaR,RandomGen);
    //if (!W.makeMove(thisWalker,deltaR, m_sqrttau))
    if (!W.makeMove(thisWalker,deltaR,SqrtTauOverMass))
    {
      H.rejectedMove(W,thisWalker);
      return;
    }
    //W.R = m_sqrttau*deltaR + thisWalker.R;
    //W.update();
    RealType logpsi(Psi.evaluateLog(W));
   // app_log()<<logpsi<< std::endl;
    RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)));
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
      H.rejectedMove(W,thisWalker);
    }
    else
    {
      RealType eloc=H.evaluate(W);
      thisWalker.R = W.R;
      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
      ++nAccept;
    }
}

//   void VMCUpdateAll::advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i)
//   {
//     int NumThreads=pclone.size();
//
//     RealType psi2_new(1.0);
//     RealType psi2_0_new(2.0*W[0]->getPropertyBase()[LOGPSI] + c_i[0]);
//     for (int ip=1; ip<NumThreads; ++ip)
//     {
//       psi2_new += std::exp(2.0*W[ip]->getPropertyBase()[LOGPSI]+ c_i[ip] -psi2_0_new);
//     }
//
//     #pragma omp parallel
//     {
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
//
//       RealType psi2_now=psi2_new;
//       RealType psi2_0_now=psi2_0_new;
//       RealType psi2_i_now=2.0*W[ip]->getPropertyBase()[LOGPSI]+ c_i[ip];
//
//     for (int iter=0; iter<nSubSteps; ++iter)
//     {
//         makeGaussRandomWithEngine(deltaR,rng_loc);
//         bool moved =  wclone[ip]->makeMove(*W[ip],deltaR, m_sqrttau);
//         if (!moved) continue;
//
//         RealType psi2_i_new(2.0*pclone[ip]->evaluateLog(*wclone[ip])+ c_i[ip]);
//
//         #pragma omp master
//        {
//          psi2_new=0.0;
//          psi2_0_new=psi2_i_new;
//        }
//
//         #pragma  omp barrier
// // #pragma flush(psi2_new,psi2_0_new)
//
//        //everybody adds the value including master
// #pragma omp critical
//        {
//          psi2_new+=std::exp(psi2_i_new-psi2_0_new);
//        }
// #pragma omp barrier
//
//        RealType prob = std::exp(psi2_0_new-psi2_0_now)*(psi2_new/psi2_now);
//         if (rng_loc() > prob) W[ip]->Age++;
//         else
//           {
//             psi2_i_now = psi2_i_new;
//             psi2_now = psi2_new;
//             psi2_0_now = psi2_0_new;
//           }
//     }
//     wclone[ip]->saveWalker(*W[ip]);
//     RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//     (*W[ip]).resetProperty(0.5*(psi2_i_now-c_i[ip]),pclone[ip]->getPhase(), eloc);
//     hclone[ip]->auxHevaluate(*wclone[ip],*W[ip]);
//     hclone[ip]->saveProperty((*W[ip]).getPropertyBase());
//     }
//   }
//
//   void VMCUpdateAll::estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//     , std::vector<MCWalkerConfiguration*>& wclone
//     , std::vector<QMCHamiltonian*>& hclone
//     , std::vector<RandomGenerator_t*>& rng
//     , std::vector<RealType>& ratio_i_0)
//   {
//     int NumThreads=pclone.size();
//
//     RealType psi2_i_now[128];
//     for(int i=0;i<NumThreads;i++)
//       psi2_i_now[i]=2.0*W[i]->getPropertyBase()[LOGPSI];
//     RealType nn(-1.0*std::log(1.0*nSubSteps));
//
//
// #pragma omp parallel
//     {
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
//
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         makeGaussRandomWithEngine(deltaR,rng_loc);
//         if (!wclone[ip]->makeMove(*W[ip],deltaR, m_sqrttau)) continue;
//         RealType psi2_i_new(2.0*pclone[ip]->evaluateLog(*wclone[ip]));
//
// #pragma omp barrier
// #pragma omp critical
//        {
//          psi2_i_now[ip]=psi2_i_new;
//        }
// #pragma omp barrier
// #pragma omp flush
// #pragma omp critical
// {
//          ratio_i_0[ip] += expl(psi2_i_now[ip]-psi2_i_now[0] + nn);
// }
// #pragma omp barrier
//     }
//
//     }
//   }

/// Constructor.
VMCUpdateAllWithDrift::VMCUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,h,rg)
{
  UpdatePbyP=false;
}

VMCUpdateAllWithDrift::~VMCUpdateAllWithDrift()
{
}

void VMCUpdateAllWithDrift::advanceWalker(Walker_t& thisWalker, bool recompute)
{
    W.loadWalker(thisWalker,false);
  //  RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    assignDrift(Tau,MassInvP,W.G,drift);
    //RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    //if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR, m_sqrttau))
    if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR,SqrtTauOverMass))
    {
      H.rejectedMove(W,thisWalker);
      return;
    }
    RealType logpsi(Psi.evaluateLog(W));
    RealType logGf = -0.5*Dot(deltaR,deltaR);
    //nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
  //  nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    assignDrift(Tau,MassInvP,W.G,drift);

    deltaR = thisWalker.R - W.R - drift;

    RealType logGb=logBackwardGF(deltaR);
    //RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);

    RealType g= std::exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
      H.rejectedMove(W,thisWalker);
    }
    else
    {
      W.saveWalker(thisWalker);
      RealType eloc=H.evaluate(W);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
      ++nAccept;
    }
}

VMCUpdateAllWithDrift::RealType VMCUpdateAllWithDrift::advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios)
{
//     Walker_t::Buffer_t& w_buffer(w1.DataSet);
  W.loadWalker(w1,false);
  std::vector<RealType> orb_ratios(4,1.0);
  int nshells(rs.size());
  std::vector<RealType> orb_phases0,orb_logs0;
  Psi.evaluateLog(W);
  RealType pt=Psi.getPhase();
  RealType lt=Psi.getLogPsi();
  Psi.getPhases(orb_phases0);
  Psi.getLogs(orb_logs0);
  int nmoved(0);
  int i(0);
  while(i<rs.size())
  {
    deltaR=0;
    while(nmoved<rs[i])
    {
      int iat=iats[nmoved];
      W.R[iat]+=dR[nmoved];
      nmoved++;
    }
    RealType logpsi(Psi.evaluateLog(W));
    std::vector<RealType> orb_phases,orb_logs;
    Psi.getPhases(orb_phases);
    Psi.getLogs(orb_logs);
    orb_ratios[0]=std::cos(pt-Psi.getPhase())*std::exp(logpsi-lt);
    orb_ratios[1]=std::cos(orb_phases[0] - orb_phases0[0])*std::exp(orb_logs[0] - orb_logs0[0]);
    orb_ratios[2]=std::cos(orb_phases[1] - orb_phases0[1])*std::exp(orb_logs[1] - orb_logs0[1]);
    orb_ratios[3]=orb_ratios[0]/orb_ratios[1];
    while(nmoved==rs[i])
    {
      ratios[i]=orb_ratios[0];
      ratios[nshells+i]=orb_ratios[1];
      ratios[nshells*2+i]=orb_ratios[2];
      ratios[nshells*3+i]=orb_ratios[3];
      i++;
      if(i==rs.size())
        break;
    }
  }
  return 0.0;
}

// void VMCUpdateAllWithDrift::advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i)
//   {
//     int NumThreads=pclone.size();
//     bool moved(false);
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
//         ParticleSet::ParticlePos_t grad_now;
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_now += std::exp(psi2_i_now[ip]-psi2_i_now[0])/psi2_now*wclone[ip]->G;
//
//         setScaledDriftPbyPandNodeCorr(m_tauovermass,grad_now,drift);
//         makeGaussRandomWithEngine(deltaR,RandomGen);
//
//         for (int ip=1; ip<NumThreads; ++ip)
//           wclone[ip]->makeMoveWithDrift(*W[ip],drift,deltaR, m_sqrttau);
//         if (!wclone[0]->makeMoveWithDrift(*W[0],drift,deltaR, m_sqrttau))
//           {
//             continue;
//           }
//
//         std::vector<RealType> psi2_i_new(NumThreads);
//         for (int ip=0; ip<NumThreads; ++ip) psi2_i_new[ip]= 2.0*pclone[ip]->evaluateLog(*wclone[ip]);
//
//         RealType psi2_new(1.0);
//         for (int ip=1; ip<NumThreads; ++ip) psi2_new += std::exp(psi2_i_new[ip]-psi2_i_new[0]);
//
//         ParticleSet::ParticlePos_t grad_new;
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_new += std::exp(psi2_i_new[ip]-psi2_i_new[0])/psi2_new*wclone[ip]->G;
//         setScaledDriftPbyPandNodeCorr(m_tauovermass,grad_new,drift);
//
//
//
//         RealType logGf = -0.5*Dot(deltaR,deltaR);
//         deltaR = W[0]->R - wclone[0]->R - drift;
//         RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
//
//         RealType p= std::exp(psi2_i_new[0]-psi2_i_now[0])*(psi2_new/psi2_now);
//         RealType g= std::exp(logGb-logGf)*p;
//
//
//         if (RandomGen() > g)
//           {
//             for (int ip=0; ip<NumThreads; ++ip) W[ip]->Age++;
//             ++nReject;
//           }
//         else
//           {
//             moved=true;
//             for (int ip=0; ip<NumThreads; ++ip) psi2_i_now[ip] = psi2_i_new[ip];
//             psi2_now = psi2_new;
//             for (int ip=0; ip<NumThreads; ++ip) wclone[ip]->saveWalker(*W[ip]);
//             ++nAccept;
//           }
//     }
// // #pragma omp parallel for
//     for (int ip=0; ip<NumThreads; ++ip)
//     {
//       Walker_t& thisWalker(*W[ip]);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       if (moved)
//         {
//           RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//           hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//           hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//         }
//     }
//     }
//


//   VMCUpdateAllSampleRN::VMCUpdateAllSampleRN(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide,
//       QMCHamiltonian& h, RandomGenerator_t& rg):
//       QMCUpdateBase(w,psi,guide,h,rg), logEpsilon(0.0)
//   {
//   }
//
//   VMCUpdateAllSampleRN::~VMCUpdateAllSampleRN()
//   {
//   }
//
//   void VMCUpdateAllSampleRN::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
//     for (; it!= it_end; ++it)
//       {
//         MCWalkerConfiguration::Walker_t& thisWalker(**it);
//         makeGaussRandomWithEngine(deltaR,RandomGen);
//         if (!W.makeMove(thisWalker,deltaR,m_sqrttau)) continue;
//         //W.R = m_sqrttau*deltaR + thisWalker.R;
//         //W.update();
//         RealType logpsi_now=thisWalker.Properties(LOGPSI);
//
//         RealType logpsi_new(Psi.evaluateLog(W));
//         RealType g= std::exp(2.0*(logpsi_new-logpsi_now))*(1+std::exp(2.0*(logEpsilon-logpsi_new)))/(1+std::exp(2.0*(logEpsilon-logpsi_now)));
//         if (RandomGen() > g)
//           {
//             thisWalker.Age++;
//             ++nReject;
//             H.rejectedMove(W,thisWalker);
//           }
//         else
//           {
//             logpsi_now = logpsi_new;
//             RealType eloc=H.evaluate(W);
//             thisWalker.R = W.R;
//             thisWalker.resetProperty(logpsi_now,Psi.getPhase(),eloc);
//             H.auxHevaluate(W,thisWalker);
//             thisWalker.Weight = 1.0/(1+std::exp(logEpsilon-logpsi_now));
//             H.saveProperty(thisWalker.getPropertyBase());
//             ++nAccept;
//           }
//       }
//   }

}


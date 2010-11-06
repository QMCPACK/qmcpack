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
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
  {

    VMCUpdateAll::VMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi,
        QMCHamiltonian& h, RandomGenerator_t& rg)
      : QMCUpdateBase(w,psi,h,rg)
      {
      }

  VMCUpdateAll::~VMCUpdateAll()
  {
  }

  void VMCUpdateAll::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
  {
    for (; it!= it_end; ++it)
      {
        MCWalkerConfiguration::Walker_t& thisWalker(**it);
        makeGaussRandomWithEngine(deltaR,RandomGen);
        if (!W.makeMove(thisWalker,deltaR, m_sqrttau)) 
        {
          H.rejectedMove(W,thisWalker);
          continue;
        }

        //W.R = m_sqrttau*deltaR + thisWalker.R;
        //W.update();

        RealType logpsi(Psi.evaluateLog(W));
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
  }

  void VMCUpdateAll::advanceCSWalkers(vector<TrialWaveFunction*>& pclone, vector<MCWalkerConfiguration*>& wclone, vector<QMCHamiltonian*>& hclone)
  {
    int NumThreads=pclone.size();
    bool moved(false);
    std::vector<RealType> psi2_i_now(NumThreads);
    RealType psi2_now=0;
    for (int ip=0; ip<NumThreads; ++ip)
    {
      psi2_i_now[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
      psi2_now += std::exp(psi2_i_now[ip]-psi2_i_now[0]);
    }
    
    for (int iter=0; iter<nSubSteps; ++iter)
    {
        makeGaussRandomWithEngine(deltaR,RandomGen);
        
        for (int ip=1; ip<NumThreads; ++ip)
          wclone[ip]->makeMove(*W[ip],deltaR, m_sqrttau);
        if (!wclone[0]->makeMove(*W[0],deltaR, m_sqrttau))
          {
            continue;
          }
          

        //W.R = m_sqrttau*deltaR + thisWalker.R;
        //W.update();
        std::vector<RealType> psi2_i_new(NumThreads);
        for (int ip=0; ip<NumThreads; ++ip) psi2_i_new[ip]= 2.0*pclone[ip]->evaluateLog(*wclone[ip]);
        RealType psi2_new(1.0);
        for (int ip=1; ip<NumThreads; ++ip) psi2_new += std::exp(psi2_i_new[ip]-psi2_i_new[0]);
        
        RealType p= std::exp(psi2_i_new[0]-psi2_i_now[0])*(psi2_new/psi2_now);

        if (RandomGen() > p)
          {
            for (int ip=0; ip<NumThreads; ++ip) W[ip]->Age++;
            ++nReject;
          }
        else
          {
            moved=true;
            for (int ip=0; ip<NumThreads; ++ip) psi2_i_now[ip] = psi2_i_new[ip];
            psi2_now = psi2_new;
            for (int ip=0; ip<NumThreads; ++ip) wclone[ip]->saveWalker(*W[ip]);
            ++nAccept;
          }
    }
// #pragma omp parallel for  
    for (int ip=0; ip<NumThreads; ++ip)
    {
      Walker_t& thisWalker(*W[ip]);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      if (moved)
        {
          RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
          thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
          hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
          hclone[ip]->saveProperty(thisWalker.getPropertyBase());
        }
    }
    }
  
  /// Constructor.
  VMCUpdateAllWithDrift::VMCUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
      QMCHamiltonian& h, RandomGenerator_t& rg):
      QMCUpdateBase(w,psi,h,rg)
  {
  }

  VMCUpdateAllWithDrift::~VMCUpdateAllWithDrift()
  {
  }

  void VMCUpdateAllWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
  {
    for (;it != it_end;++it)
      {
        MCWalkerConfiguration::Walker_t& thisWalker(**it);
        W.loadWalker(thisWalker,false);
        RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
        
        makeGaussRandomWithEngine(deltaR,RandomGen);
        if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR, m_sqrttau))
        {
          H.rejectedMove(W,thisWalker);
          continue;
        }
        //
        //W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
        //W.update();
        RealType logpsi(Psi.evaluateLog(W));
        RealType logGf = -0.5*Dot(deltaR,deltaR);
        // setScaledDrift(m_tauovermass,W.G,drift);
        nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);

        //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
        deltaR = thisWalker.R - W.R - drift;
        RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);

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
  }

void VMCUpdateAllWithDrift::advanceCSWalkers(vector<TrialWaveFunction*>& pclone, vector<MCWalkerConfiguration*>& wclone, vector<QMCHamiltonian*>& hclone)
  {
    int NumThreads=pclone.size();
    bool moved(false);
    std::vector<RealType> psi2_i_now(NumThreads);
    RealType psi2_now=0;
    for (int ip=0; ip<NumThreads; ++ip)
    {
      psi2_i_now[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
      psi2_now += std::exp(psi2_i_now[ip]-psi2_i_now[0]);
    }
    
    for (int iter=0; iter<nSubSteps; ++iter)
    {
        ParticleSet::ParticlePos_t grad_now;
        for (int ip=0; ip<NumThreads; ++ip)
          grad_now += std::exp(psi2_i_now[ip]-psi2_i_now[0])/psi2_now*wclone[ip]->G;
        
        setScaledDriftPbyPandNodeCorr(m_tauovermass,grad_now,drift);
        makeGaussRandomWithEngine(deltaR,RandomGen);
        
        for (int ip=1; ip<NumThreads; ++ip)
          wclone[ip]->makeMoveWithDrift(*W[ip],drift,deltaR, m_sqrttau);
        if (!wclone[0]->makeMoveWithDrift(*W[0],drift,deltaR, m_sqrttau))
          {
            continue;
          }
        
        std::vector<RealType> psi2_i_new(NumThreads);
        for (int ip=0; ip<NumThreads; ++ip) psi2_i_new[ip]= 2.0*pclone[ip]->evaluateLog(*wclone[ip]);
        
        RealType psi2_new(1.0);
        for (int ip=1; ip<NumThreads; ++ip) psi2_new += std::exp(psi2_i_new[ip]-psi2_i_new[0]);
        
        ParticleSet::ParticlePos_t grad_new;
        for (int ip=0; ip<NumThreads; ++ip)
          grad_new += std::exp(psi2_i_new[ip]-psi2_i_new[0])/psi2_new*wclone[ip]->G;
        setScaledDriftPbyPandNodeCorr(m_tauovermass,grad_new,drift);
        
        
        
        RealType logGf = -0.5*Dot(deltaR,deltaR);
        deltaR = W[0]->R - wclone[0]->R - drift;
        RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);

        RealType p= std::exp(psi2_i_new[0]-psi2_i_now[0])*(psi2_new/psi2_now);
        RealType g= std::exp(logGb-logGf)*p; 
        

        if (RandomGen() > g)
          {
            for (int ip=0; ip<NumThreads; ++ip) W[ip]->Age++;
            ++nReject;
          }
        else
          {
            moved=true;
            for (int ip=0; ip<NumThreads; ++ip) psi2_i_now[ip] = psi2_i_new[ip];
            psi2_now = psi2_new;
            for (int ip=0; ip<NumThreads; ++ip) wclone[ip]->saveWalker(*W[ip]);
            ++nAccept;
          }
    }
// #pragma omp parallel for  
    for (int ip=0; ip<NumThreads; ++ip)
    {
      Walker_t& thisWalker(*W[ip]);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      if (moved)
        {
          RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
          thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
          hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
          hclone[ip]->saveProperty(thisWalker.getPropertyBase());
        }
    }
    }
    


  VMCUpdateAllSampleRN::VMCUpdateAllSampleRN(MCWalkerConfiguration& w, TrialWaveFunction& psi,
      QMCHamiltonian& h, RandomGenerator_t& rg):
      QMCUpdateBase(w,psi,h,rg), logEpsilon(0.0)
  {
  }

  VMCUpdateAllSampleRN::~VMCUpdateAllSampleRN()
  {
  }

  void VMCUpdateAllSampleRN::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
  {
    for (; it!= it_end; ++it)
      {
        MCWalkerConfiguration::Walker_t& thisWalker(**it);
        makeGaussRandomWithEngine(deltaR,RandomGen);
        if (!W.makeMove(thisWalker,deltaR,m_sqrttau)) continue;
        //W.R = m_sqrttau*deltaR + thisWalker.R;
        //W.update();
        RealType logpsi_now=thisWalker.Properties(LOGPSI);

        RealType logpsi_new(Psi.evaluateLog(W));
        RealType g= std::exp(2.0*(logpsi_new-logpsi_now))*(1+std::exp(2.0*(logEpsilon-logpsi_new)))/(1+std::exp(2.0*(logEpsilon-logpsi_now)));
        if (RandomGen() > g)
          {
            thisWalker.Age++;
            ++nReject;
            H.rejectedMove(W,thisWalker);
          }
        else
          {
            logpsi_now = logpsi_new;
            RealType eloc=H.evaluate(W);
            thisWalker.R = W.R;
            thisWalker.resetProperty(logpsi_now,Psi.getPhase(),eloc);
            H.auxHevaluate(W,thisWalker);
            thisWalker.Weight = 1.0/(1+std::exp(logEpsilon-logpsi_now));
            H.saveProperty(thisWalker.getPropertyBase());
            ++nAccept;
          }
      }
  }

}

/***************************************************************************
 * $RCSfile: VMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdateAll.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/

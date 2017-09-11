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
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/CorrelatedSampling/CSVMCUpdateAll.h"
//#include "Utilities/OhmmsInfo.h"
//#include "Particle/MCWalkerConfiguration.h"
//#include "Particle/HDFWalkerIO.h"
//#include "ParticleBase/ParticleUtility.h"
//#include "ParticleBase/RandomSeqGenerator.h"
//#include "ParticleBase/ParticleAttribOps.h"
//#include "Message/Communicate.h"
//#include "Estimators/MultipleEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

/// Constructor.
CSVMCUpdateAll::CSVMCUpdateAll(MCWalkerConfiguration& w,
                                std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& h, RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{ 
  UpdatePbyP=false;	
}

void CSVMCUpdateAll::advanceWalker(Walker_t& thisWalker, bool recompute)
{
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
   
    RealType tau_over_mass = std::sqrt(Tau*MassInvS[0]);
    
    if (!W.makeMove(thisWalker,deltaR,tau_over_mass))
    {
      for (int ipsi=0; ipsi<nPsi; ipsi++) H1[ipsi]->rejectedMove(W,thisWalker);
      return;
    }
    
    //Evaluate Psi and graidients and laplacians
    //\f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      logpsi[ipsi]=Psi1[ipsi]->evaluateLog(W);
      Psi1[ipsi]->L=W.L;
      Psi1[ipsi]->G=W.G;

     *G1[ipsi]=W.G;
     *L1[ipsi]=W.L;    
      sumratio[ipsi]=1.0;
    }
    // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
    for(int ipsi=0; ipsi< nPsi-1; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType ratioij=avgNorm[ipsi]/avgNorm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
        sumratio[ipsi] += ratioij;
        sumratio[jpsi] += 1.0/ratioij;
      }
      
    }
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      invsumratio[ipsi]=1.0/sumratio[ipsi];
      cumNorm[ipsi]+=invsumratio[ipsi];
    }
    
       for(int ipsi=0; ipsi<nPsi; ipsi++) cumNorm[ipsi]+=invsumratio[ipsi];
       
    RealType g = sumratio[0]/thisWalker.Multiplicity*
                 std::exp(2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
    
    if(Random() > g)
    {
      thisWalker.Age++;
      ++nReject;
      for (int ipsi=0; ipsi<nPsi; ipsi++) H1[ipsi]->rejectedMove(W,thisWalker);
    }
    else
    {
      thisWalker.Age=0;
      thisWalker.Multiplicity=sumratio[0];
      thisWalker.R = W.R;
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        W.L=*L1[ipsi];
        W.G=*G1[ipsi];

        RealType et = H1[ipsi]->evaluate(W);
        thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
        thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
        thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
        thisWalker.Properties(ipsi,LOCALENERGY)=et;
         H1[ipsi]->auxHevaluate(W,thisWalker);
         H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      }
      ++nAccept;
    }

}


CSVMCUpdateAllWithDrift::CSVMCUpdateAllWithDrift(MCWalkerConfiguration& w,
                                std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& h, RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{ 
  UpdatePbyP=false;
}

void CSVMCUpdateAllWithDrift::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  //create a 3N-Dimensional Gaussian with variance=1
  W.loadWalker(thisWalker,false);
  assignDrift(Tau,MassInvP,W.G,drift);
  makeGaussRandomWithEngine(deltaR,RandomGen);

  Walker_t::ParticleGradient_t cumGrad(W.G);
  cumGrad=0.0;

  RealType tau_over_mass = std::sqrt(Tau*MassInvS[0]);

  if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR,SqrtTauOverMass))
  {
    for (int ipsi=1; ipsi<nPsi; ipsi++) H1[ipsi]->rejectedMove(W,thisWalker);
    return;
  }



  RealType logGf = -0.5*Dot(deltaR,deltaR);

  //Evaluate Psi and graidients and laplacians
  //\f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    logpsi[ipsi]=Psi1[ipsi]->evaluateLog(W);
    Psi1[ipsi]->L=W.L;
    Psi1[ipsi]->G=W.G;
    *G1[ipsi]=W.G;
    *L1[ipsi]=W.L;    
    sumratio[ipsi]=1.0;
  }
  // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
  for(int ipsi=0; ipsi< nPsi-1; ipsi++)
  {
    for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
    {
      RealType ratioij=avgNorm[ipsi]/avgNorm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
      sumratio[ipsi] += ratioij;
      sumratio[jpsi] += 1.0/ratioij;
    }

  }
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    invsumratio[ipsi]=1.0/sumratio[ipsi];
    cumGrad+=Psi1[ipsi]->G*static_cast<Walker_t::ParticleValue_t>(invsumratio[ipsi]);
  }

    for(int ipsi=0; ipsi<nPsi; ipsi++) cumNorm[ipsi]+=invsumratio[ipsi];

  assignDrift(Tau,MassInvP,cumGrad,drift);

  deltaR = thisWalker.R - W.R - drift;

  RealType logGb=logBackwardGF(deltaR);

  RealType g = sumratio[0]/thisWalker.Multiplicity*
    std::exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));

  if(Random() > g)
  {
    thisWalker.Age++;
    ++nReject;
    for (int ipsi=1; ipsi<nPsi; ipsi++) H1[ipsi]->rejectedMove(W,thisWalker);
  }
  else
  {
    thisWalker.Age=0;
    thisWalker.Multiplicity=sumratio[0];
    thisWalker.R = W.R;
    thisWalker.G=cumGrad;
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      W.L=*L1[ipsi];
      W.G=*G1[ipsi];
      RealType et = H1[ipsi]->evaluate(W);
      thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
      thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
      thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
      thisWalker.Properties(ipsi,LOCALENERGY)=et;
      H1[ipsi]->auxHevaluate(W,thisWalker);
      H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
    }
    ++nAccept;
  }
}

}


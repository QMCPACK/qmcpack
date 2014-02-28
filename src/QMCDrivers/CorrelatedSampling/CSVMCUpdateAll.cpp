//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Simone Chiesa
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
                                vector<TrialWaveFunction*>& psi, vector<QMCHamiltonian*>& h, RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{ 
  UpdatePbyP=false;	
}

/**  Advance all the walkers one timstep.
 */
void CSVMCUpdateAll::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  
  for (; it != it_end; ++it)
  {
    MCWalkerConfiguration::Walker_t &thisWalker(**it);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
   
    RealType tau_over_mass = std::sqrt(Tau*MassInvS[0]);
   // app_log()<<"tau_over_mass= "<<tau_over_mass<<endl;
    //app_log()<<"deltaR = "<<deltaR<<endl;
    //if (!W.makeMove(thisWalker,deltaR, m_sqrttau))
    
    if (!W.makeMove(thisWalker,deltaR,tau_over_mass))
    {
      for (int ipsi=0; ipsi<nPsi; ipsi++) H1[ipsi]->rejectedMove(W,thisWalker);
      continue;
    }
    
   // if(useDrift)
  //  {
  ///    //forward green function
  ///    RealType logGf = -0.5*Dot(deltaR,deltaR);
  ///    PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
  ///    for(int ipsi=1; ipsi< nPsi ; ipsi++)
  ///    {
  ///      PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
  ///    }
  ///    setScaledDrift(Tau,drift);
  ///    //backward green function
  ///   deltaR = thisWalker.R - W.R - drift;
  ///    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
  ///    g *= std::exp(logGb-logGf);
  ///  }
    
  ///  if(useDrift)
   ///   W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
  ///  else
  ///    W.R = m_sqrttau*deltaR + thisWalker.R;
  
  
    //update the distance table associated with W
    //DistanceTable::update(W);
  ////  W.update();
    //Evaluate Psi and graidients and laplacians
    //\f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      //W.L=0;
      //W.G=0;
     // app_log()<<"psi0 "<<Psi.evaluateLog(W)<<endl;
      logpsi[ipsi]=Psi1[ipsi]->evaluateLog(W);
      Psi1[ipsi]->L=W.L;
      Psi1[ipsi]->G=W.G;
   //   W.L=0;
    //  W.G=0;
     *G1[ipsi]=W.G;
     *L1[ipsi]=W.L;    
//      app_log()<<"Psiname="<<Psi1[ipsi]->getName()<<" H name = "<<H1[ipsi]->getName()<<endl;
    //  app_log()<<ipsi<<" logpsi "<<logpsi[ipsi]<<endl;
    //  app_log()<<ipsi<<" L "<<W.L<<endl;
   //   app_log()<<ipsi<<" Lsave"<<*L1[ipsi]<<endl;
   //   app_log()<<ipsi<<" G "<<W.G<<endl;
    //  app_log()<<ipsi<<" Gsave "<<*G1[ipsi]<<endl;
      sumratio[ipsi]=1.0;
    }
  //  app_log()<<endl;
    // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
    for(int ipsi=0; ipsi< nPsi-1; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType ratioij=avgNorm[ipsi]/avgNorm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
       // app_log()<<ipsi<<" "<<jpsi<<" "<<ratioij<<" "<<logpsi[ipsi]<<" "<<logpsi[jpsi]<<" "<<avgNorm[ipsi]<<" "<<avgNorm[jpsi]<<endl;
        sumratio[ipsi] += ratioij;
        sumratio[jpsi] += 1.0/ratioij;
      }
      
    }
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      invsumratio[ipsi]=1.0/sumratio[ipsi];
      cumNorm[ipsi]+=invsumratio[ipsi];
    }
    
    if (measure==true)
       for(int ipsi=0; ipsi<nPsi; ipsi++) cumNorm[ipsi]+=invsumratio[ipsi];
       
    RealType g = sumratio[0]/thisWalker.Multiplicity*
                 std::exp(2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
 ///   if(useDrift)
 ///   {
  ///    //forward green function
  ///    RealType logGf = -0.5*Dot(deltaR,deltaR);
  ///    PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
  ///    for(int ipsi=1; ipsi< nPsi ; ipsi++)
  ///    {
  ///      PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
  ///    }
  ///    setScaledDrift(Tau,drift);
  ///    //backward green function
  ///   deltaR = thisWalker.R - W.R - drift;
  ///    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
  ///    g *= std::exp(logGb-logGf);
  ///  }
    //Original
    //RealType g = Properties(SUMRATIO)/thisWalker.Properties(SUMRATIO)*
    //	exp(logGb-logGf+2.0*(Properties(LOGPSI)-thisWalker.Properties(LOGPSI)));
    //Reuse Multiplicity to store the sumratio[0]
    //This is broken up into two pieces
    //RealType g = sumratio[0]/thisWalker.Multiplicity*
    // 	std::exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
    
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
     // thisWalker.Drift = drift;
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
      //  W.L=Psi1[ipsi]->L;
       // W.G=Psi1[ipsi]->G;
       // app_log()<<ipsi<<" L "<<W.L<<endl;
       // app_log()<<ipsi<<" G "<<W.G<<endl;
        W.L=*L1[ipsi];
        W.G=*G1[ipsi];
     //   W.L=0;
     //   W.G=0;
        RealType et = H1[ipsi]->evaluate(W);
        thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
        thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
        thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
       // cumNorm[ipsi]+=invsumratio[ipsi];
        thisWalker.Properties(ipsi,LOCALENERGY)=et;
        //multiEstimator->updateSample(iwlk,ipsi,et,invsumratio[ipsi]);
        //app_log()<<thisWalker.Properties(ipsi,UMBRELLAWEIGHT)<<" "<<endl;
        H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      }
    //  updateNorms();
       // app_log()<<endl;
      ++nAccept;
    }

 // }
 }
}


CSVMCUpdateAllWithDrift::CSVMCUpdateAllWithDrift(MCWalkerConfiguration& w,
                                vector<TrialWaveFunction*>& psi, vector<QMCHamiltonian*>& h, RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{ 
  UpdatePbyP=false;
}

void CSVMCUpdateAllWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  for (; it != it_end; ++it)
  {
    MCWalkerConfiguration::Walker_t &thisWalker(**it);
    //create a 3N-Dimensional Gaussian with variance=1
    W.loadWalker(thisWalker,false);
    assignDrift(Tau,MassInvP,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    
    Walker_t::ParticleGradient_t cumGrad(W.G);
    cumGrad=0.0;
    
    RealType tau_over_mass = std::sqrt(Tau*MassInvS[0]);
   // app_log()<<"tau_over_mass= "<<tau_over_mass<<endl;
    //app_log()<<"deltaR = "<<deltaR<<endl;
    //if (!W.makeMove(thisWalker,deltaR, m_sqrttau))
    
    if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR,SqrtTauOverMass))
    {
      for (int ipsi=1; ipsi<nPsi; ipsi++) H1[ipsi]->rejectedMove(W,thisWalker);
      // H.rejectedMove(W,thisWalker);
      continue;
    }
    
    
    
    RealType logGf = -0.5*Dot(deltaR,deltaR);
  
    
   // if(useDrift)
  //  {
  ///    //forward green function
  ///    RealType logGf = -0.5*Dot(deltaR,deltaR);
  ///    PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
  ///    for(int ipsi=1; ipsi< nPsi ; ipsi++)
  ///    {
  ///      PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
  ///    }
  ///    setScaledDrift(Tau,drift);
  ///    //backward green function
  ///   deltaR = thisWalker.R - W.R - drift;
  ///    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
  ///    g *= std::exp(logGb-logGf);
  ///  }
    
  ///  if(useDrift)
   ///   W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
  ///  else
  ///    W.R = m_sqrttau*deltaR + thisWalker.R;
  
  
    //update the distance table associated with W
    //DistanceTable::update(W);
  ////  W.update();
    //Evaluate Psi and graidients and laplacians
    //\f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      //W.L=0;
      //W.G=0;
     // app_log()<<"psi0 "<<Psi.evaluateLog(W)<<endl;
      logpsi[ipsi]=Psi1[ipsi]->evaluateLog(W);
      Psi1[ipsi]->L=W.L;
      Psi1[ipsi]->G=W.G;
   //   W.L=0;
    //  W.G=0;
     *G1[ipsi]=W.G;
     *L1[ipsi]=W.L;    
//      app_log()<<"Psiname="<<Psi1[ipsi]->getName()<<" H name = "<<H1[ipsi]->getName()<<endl;
    //  app_log()<<ipsi<<" logpsi "<<logpsi[ipsi]<<endl;
    //  app_log()<<ipsi<<" L "<<W.L<<endl;
   //   app_log()<<ipsi<<" Lsave"<<*L1[ipsi]<<endl;
   //   app_log()<<ipsi<<" G "<<W.G<<endl;
    //  app_log()<<ipsi<<" Gsave "<<*G1[ipsi]<<endl;
      sumratio[ipsi]=1.0;
    }
  //  app_log()<<endl;
    // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
    for(int ipsi=0; ipsi< nPsi-1; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType ratioij=avgNorm[ipsi]/avgNorm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
       // app_log()<<ipsi<<" "<<jpsi<<" "<<ratioij<<" "<<logpsi[ipsi]<<" "<<logpsi[jpsi]<<" "<<avgNorm[ipsi]<<" "<<avgNorm[jpsi]<<endl;
        sumratio[ipsi] += ratioij;
        sumratio[jpsi] += 1.0/ratioij;
      }
      
    }
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      invsumratio[ipsi]=1.0/sumratio[ipsi];
      cumGrad+=Psi1[ipsi]->G*invsumratio[ipsi];
    }
    
    if (measure==true)
       for(int ipsi=0; ipsi<nPsi; ipsi++) cumNorm[ipsi]+=invsumratio[ipsi];
       
    assignDrift(Tau,MassInvP,W.G,drift);

    deltaR = thisWalker.R - W.R - drift;

    RealType logGb=logBackwardGF(deltaR);
       
    RealType g = sumratio[0]/thisWalker.Multiplicity*
                 std::exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
 ///   if(useDrift)
 ///   {
  ///    //forward green function
  ///    RealType logGf = -0.5*Dot(deltaR,deltaR);
  ///    PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
  ///    for(int ipsi=1; ipsi< nPsi ; ipsi++)
  ///    {
  ///      PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
  ///    }
  ///    setScaledDrift(Tau,drift);
  ///    //backward green function
  ///   deltaR = thisWalker.R - W.R - drift;
  ///    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
  ///    g *= std::exp(logGb-logGf);
  ///  }
    //Original
    //RealType g = Properties(SUMRATIO)/thisWalker.Properties(SUMRATIO)*
    //	exp(logGb-logGf+2.0*(Properties(LOGPSI)-thisWalker.Properties(LOGPSI)));
    //Reuse Multiplicity to store the sumratio[0]
    //This is broken up into two pieces
    //RealType g = sumratio[0]/thisWalker.Multiplicity*
    // 	std::exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
    
    if(Random() > g)
    {
      thisWalker.Age++;
      ++nReject;
    //  H.rejectedMove(W,thisWalker);
      for (int ipsi=1; ipsi<nPsi; ipsi++) H1[ipsi]->rejectedMove(W,thisWalker);
    }
    else
    {
      thisWalker.Age=0;
      thisWalker.Multiplicity=sumratio[0];
      thisWalker.R = W.R;
      thisWalker.G=cumGrad;
     // thisWalker.Drift = drift;
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
      //  W.L=Psi1[ipsi]->L;
       // W.G=Psi1[ipsi]->G;
       // app_log()<<ipsi<<" L "<<W.L<<endl;
       // app_log()<<ipsi<<" G "<<W.G<<endl;
        W.L=*L1[ipsi];
        W.G=*G1[ipsi];
     //   W.L=0;
     //   W.G=0;
        RealType et = H1[ipsi]->evaluate(W);
        thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
        thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
        thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
       // cumNorm[ipsi]+=invsumratio[ipsi];
        thisWalker.Properties(ipsi,LOCALENERGY)=et;
        //multiEstimator->updateSample(iwlk,ipsi,et,invsumratio[ipsi]);
        //app_log()<<thisWalker.Properties(ipsi,UMBRELLAWEIGHT)<<" "<<endl;
        H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      }
    //  updateNorms();
       // app_log()<<endl;
      ++nAccept;
    }
 }
}

}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1593 $   $Date: 2007-01-04 17:23:27 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultiple.cpp 1593 2007-01-04 23:23:27Z jnkim $
 ***************************************************************************/

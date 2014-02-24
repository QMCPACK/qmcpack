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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/CorrelatedSampling/CSUpdateBase.h"
#include "Estimators/CSEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"
#include "Utilities/IteratorUtility.h"
#include <numeric>

namespace qmcplusplus
{
CSUpdateBase::CSUpdateBase(MCWalkerConfiguration& w,
                           vector<TrialWaveFunction*>& psipool, vector<QMCHamiltonian*>& hpool, RandomGenerator_t& rg):
  QMCUpdateBase(w,*psipool[0],*hpool[0],rg), nPsi(0), useDriftOption("no"), H1(hpool), Psi1(psipool)
{
  myParams.add(useDriftOption,"useDrift","string");
  
}

CSUpdateBase::~CSUpdateBase()
{
  delete_iter(G1.begin(),G1.end());
  delete_iter(L1.begin(),L1.end());
}

void CSUpdateBase::resizeWorkSpace(int nw,int nptcls)
{
  if(logpsi.size())
    return;
  logpsi.resize(nPsi);
  sumratio.resize(nPsi);
  invsumratio.resize(nPsi);
  avgNorm.resize(nPsi,1.0);
  logNorm.resize(nPsi,0.0);
  cumNorm.resize(nPsi,0.0);
  avgWeight.resize(nPsi,1.0);
  instRij.resize(nPsi*(nPsi-1)/2);
  ratioIJ.resize(nw,nPsi*(nPsi-1)/2);
  dG.resize(nptcls);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    Psi1[ipsi]->G.resize(nptcls);
    Psi1[ipsi]->L.resize(nptcls);
    G1.push_back(new ParticleSet::ParticleGradient_t(nptcls));
    L1.push_back(new ParticleSet::ParticleLaplacian_t(nptcls));
  }
}

void CSUpdateBase::updateNorms()
{
 //for(int ipsi=0; ipsi< nPsi; ipsi++)
 //   cumNorm[ipsi]+=multiEstimator->getUmbrellaWeight(ipsi);
  //if(block==(equilBlocks-1) || block==(nBlocks-1)){
//	  app_log()<<"Inside UpdateNorm\n";
  RealType winv=1.0/double(std::accumulate(cumNorm.begin(), cumNorm.end(),0.0));
  for(int ipsi=0; ipsi< nPsi; ipsi++)
  {
    avgNorm[ipsi]=cumNorm[ipsi]*winv;
   // avgNorm[ipsi]=0.5;
    logNorm[ipsi]=std::log(avgNorm[ipsi]);
   // app_log()<<ipsi<<" "<<avgNorm[ipsi]<<" "<<logNorm[ipsi]<<" "<<winv<<endl;
    cumNorm[ipsi]=0;
  }
  
  //}
}
void CSUpdateBase::updateAvgWeights()
{
	  RealType winv=1.0/double(std::accumulate(cumNorm.begin(), cumNorm.end(),0.0));
  for(int ipsi=0; ipsi< nPsi; ipsi++)
  {
    avgWeight[ipsi]=cumNorm[ipsi]*winv;
   //  app_log()<<ipsi<<" "<<avgWeight[ipsi]<<endl;
   // avgNorm[ipsi]=0.5;
    cumNorm[ipsi]=0;
  }
}

void CSUpdateBase::initCSWalkers(WalkerIter_t it, WalkerIter_t it_end,
                                 bool resetNorms)
{
  nPsi=Psi1.size();
  useDrift=(useDriftOption=="yes");
  if(nPsi ==0)
  {
    app_error() << "  CSUpdateBase::initCSWalkers fails. Empyty Psi/H pairs" << endl;
    abort();//FIX_ABORT
  }
  int nw = it_end-it;//W.getActiveWalkers();
  resizeWorkSpace(nw,W.getTotalNum());
  if(resetNorms)
    logNorm.resize(nPsi,0.0);
  for(int ipsi=0; ipsi< nPsi; ipsi++)
    avgNorm[ipsi]=std::exp(logNorm[ipsi]);
  int iw(0);
  while(it != it_end)
  {
    Walker_t& thisWalker(**it);
    W.R = thisWalker.R;
    
    W.update();
    //evalaute the wavefunction and hamiltonian
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      logpsi[ipsi]=Psi1[ipsi]->evaluateLog(W);
      Psi1[ipsi]->G=W.G;
      thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
      RealType e=thisWalker.Properties(ipsi,LOCALENERGY)=H1[ipsi]->evaluate(W);
      H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      sumratio[ipsi]=1.0;
    }
    //Check SIMONE's note
    //Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
    int indexij(0);
    RealType* restrict rPtr=ratioIJ[iw];
    for(int ipsi=0; ipsi< nPsi-1; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType r=std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
        rPtr[indexij++]=r*avgNorm[ipsi]/avgNorm[jpsi];
        sumratio[ipsi] += r;
        sumratio[jpsi] += 1.0/r;
      }
    }
    //Re-use Multiplicity as the sumratio
    thisWalker.Multiplicity=sumratio[0];
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      thisWalker.Properties(ipsi,UMBRELLAWEIGHT)
      = invsumratio[ipsi] =1.0/sumratio[ipsi];
      cumNorm[ipsi]+=1.0/sumratio[ipsi];
    }
    //DON't forget DRIFT!!!
 ///   thisWalker.Drift=0.0;
 ///   if(useDrift)
  ///  {
 ///     for(int ipsi=0; ipsi< nPsi; ipsi++)
 ///       PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,thisWalker.Drift);
 ///     setScaledDrift(Tau,thisWalker.Drift);
  ///  }
    ++it;
    ++iw;
  }
}

void CSUpdateBase::initCSWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end,
                                        bool resetNorms)
{
  nPsi=Psi1.size();
  useDrift=(useDriftOption=="yes");
  if(nPsi ==0)
  {
    app_error() << "  CSUpdateBase::initCSWalkers fails. Empyty Psi/H pairs" << endl;
    abort();//FIX_ABORT
  }
  int nw = it_end-it;//W.getActiveWalkers();
  resizeWorkSpace(nw,W.getTotalNum());
  if(resetNorms)
    logNorm.resize(nPsi,0.0);
  for(int ipsi=0; ipsi< nPsi; ipsi++)
    avgNorm[ipsi]=std::exp(logNorm[ipsi]);
  int iw=0;
  while(it != it_end)
  {
    Walker_t& thisWalker(**it);
    thisWalker.DataSet.clear();
    thisWalker.DataSet.rewind();
   ///W.registerData(thisWalker,(*it)->DataSet);
    //evalaute the wavefunction and hamiltonian
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      //Need to modify the return value of OrbitalBase::registerData
      logpsi[ipsi]=Psi1[ipsi]->registerData(W,(*it)->DataSet);
      Psi1[ipsi]->G=W.G;
      thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
      thisWalker.Properties(ipsi,LOCALENERGY)=H1[ipsi]->evaluate(W);
      H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      sumratio[ipsi]=1.0;
    }
    //Check SIMONE's note
    //Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
    int indexij(0);
    RealType *rPtr=ratioIJ[iw];
    for(int ipsi=0; ipsi< nPsi-1; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++, indexij++)
      {
        RealType r=std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
        //rPtr[indexij++]=r*avgNorm[ipsi]/avgNorm[jpsi];
        rPtr[indexij]=r*avgNorm[ipsi]/avgNorm[jpsi];
        sumratio[ipsi] += r;
        sumratio[jpsi] += 1.0/r;
      }
    }
    //Re-use Multiplicity as the sumratio
    thisWalker.Multiplicity=sumratio[0];
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      thisWalker.Properties(ipsi,UMBRELLAWEIGHT)
      = invsumratio[ipsi] =1.0/sumratio[ipsi];
    }
    //DON't forget DRIFT!!!
///    thisWalker.Drift=0.0;
///    if(useDrift)
///    {
///      for(int ipsi=0; ipsi< nPsi; ipsi++)
///        PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,thisWalker.Drift);
///      setScaledDrift(Tau,thisWalker.Drift);
///    }
    ++it;
    ++iw;
  }
}

void CSUpdateBase::updateCSWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  int iw=0;
  while(it != it_end)
  {
    Walker_t& thisWalker(**it);
    Walker_t::Buffer_t& w_buffer((*it)->DataSet);
   //app_log()<<"DAMN.  YOU FOUND ME.  (updateCSWalkers called)\n";
   // w_buffer.rewind();
  //  W.updateBuffer(**it,w_buffer);
    //evalaute the wavefunction and hamiltonian
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      //Need to modify the return value of OrbitalBase::registerData
      logpsi[ipsi]=Psi1[ipsi]->updateBuffer(W,(*it)->DataSet);
      Psi1[ipsi]->G=W.G;
      thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
      thisWalker.Properties(ipsi,LOCALENERGY)=H1[ipsi]->evaluate(W);
      H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      sumratio[ipsi]=1.0;
    }
    int indexij(0);
    RealType *rPtr=ratioIJ[iw];
    for(int ipsi=0; ipsi< nPsi-1; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType r=std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
        rPtr[indexij++]=r*avgNorm[ipsi]/avgNorm[jpsi];
        sumratio[ipsi] += r;
        sumratio[jpsi] += 1.0/r;
      }
    }
    //Re-use Multiplicity as the sumratio
    thisWalker.Multiplicity=sumratio[0];
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      thisWalker.Properties(ipsi,UMBRELLAWEIGHT)
      = invsumratio[ipsi] =1.0/sumratio[ipsi];
    }
    //DON't forget DRIFT!!!
///    thisWalker.Drift=0.0;
///    if(useDrift)
///    {
///      for(int ipsi=0; ipsi< nPsi; ipsi++)
///        PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,thisWalker.Drift);
///      setScaledDrift(Tau,thisWalker.Drift);
///    }
    ++it;
    ++iw;
  }
}


}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1593 $   $Date: 2007-01-04 17:23:27 -0600 (Thu, 04 Jan 2007) $
 * $Id: VMCMultiple.cpp 1593 2007-01-04 23:23:27Z jnkim $
 ***************************************************************************/

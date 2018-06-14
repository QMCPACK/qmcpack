//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "Platforms/sysutil.h"
#include "QMCDrivers/CorrelatedSampling/CSUpdateBase.h"
#include "Estimators/CSEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"
#include "Utilities/IteratorUtility.h"
#include <numeric>

namespace qmcplusplus
{
CSUpdateBase::CSUpdateBase(MCWalkerConfiguration& w,
                           std::vector<TrialWaveFunction*>& psipool, std::vector<QMCHamiltonian*>& hpool, RandomGenerator_t& rg):
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
  ratio.resize(nPsi);
  sumratio.resize(nPsi);
  invsumratio.resize(nPsi);
  avgNorm.resize(nPsi,1.0);
  logNorm.resize(nPsi,0.0);
  cumNorm.resize(nPsi,0.0);
  avgWeight.resize(nPsi,1.0);
  instRij.resize(nPsi*(nPsi-1)/2);
  ratioIJ.resize(nw,nPsi*(nPsi-1)/2);
  RatioIJ.resize(nPsi,nPsi);
  dG.resize(nptcls);

  g1_new.resize(nPsi);
  g1_old.resize(nPsi);

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
  RealType winv=1.0/double( std::accumulate(cumNorm.begin(), cumNorm.end(),0.0));
  for(int ipsi=0; ipsi< nPsi; ipsi++)
  {
    avgNorm[ipsi]=cumNorm[ipsi]*winv;
   // avgNorm[ipsi]=0.5;
    logNorm[ipsi]=std::log(avgNorm[ipsi]);
   // app_log()<<ipsi<<" "<<avgNorm[ipsi]<<" "<<logNorm[ipsi]<<" "<<winv<< std::endl;
    cumNorm[ipsi]=0;
  }
  
  //}
}
void CSUpdateBase::updateAvgWeights()
{
	  RealType winv=1.0/double( std::accumulate(cumNorm.begin(), cumNorm.end(),0.0));
  for(int ipsi=0; ipsi< nPsi; ipsi++)
  {
    avgWeight[ipsi]=cumNorm[ipsi]*winv;
   //  app_log()<<ipsi<<" "<<avgWeight[ipsi]<< std::endl;
   // avgNorm[ipsi]=0.5;
    cumNorm[ipsi]=0;
  }
}

void CSUpdateBase::computeSumRatio(const std::vector<RealType> & lpsi,
                                   const std::vector<RealType> & avnorm,
                                         std::vector<RealType> & sumrat)
{
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    sumrat[ipsi]=1.0;

  for(int ipsi=0; ipsi< nPsi-1; ipsi++)
  {
    for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType ratioij=avnorm[ipsi]/avnorm[jpsi]*std::exp(2.0*(lpsi[jpsi]-lpsi[ipsi]));
        sumrat[ipsi] += ratioij;
        sumrat[jpsi] += 1.0/ratioij;
      }
  }
}

//This function assumes logpsi and avgnorm are available.  it returns
//ratioij and sumratio.
void CSUpdateBase::computeSumRatio(const std::vector<RealType> & lpsi,
                                   const std::vector<RealType> & avnorm,
                                         Matrix<RealType> & Ratioij,
                                         std::vector<RealType> & sumrat)
{
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    sumrat[ipsi]=1.0;
    Ratioij[ipsi][ipsi]=1.0;
  }

  for(int ipsi=0; ipsi< nPsi-1; ipsi++)
  {
    for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        Ratioij[ipsi][jpsi]=avnorm[ipsi]/avnorm[jpsi]*std::exp(2.0*(lpsi[jpsi]-lpsi[ipsi]));
        Ratioij[jpsi][ipsi]=1.0/Ratioij[ipsi][jpsi];
        sumrat[ipsi] += Ratioij[ipsi][jpsi];
        sumrat[jpsi] += Ratioij[jpsi][ipsi];
      }
  }
}

//This function assumes that Ratioij has been given, then computes sumratio.
void CSUpdateBase::computeSumRatio(const Matrix<RealType> & Ratioij,
                                    std::vector<RealType> & sumrat)
{
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    sumrat[ipsi]=1.0;

  for(int ipsi=0; ipsi< nPsi-1; ipsi++)
  {
    for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        sumrat[ipsi] += Ratioij[ipsi][jpsi];
        sumrat[jpsi] += Ratioij[jpsi][ipsi];
      }
  }
}

//this function takes a matrix containing |psi_i(r)/psi_j(r)|^2, and a vector
//containing |\psi_i(r')/\psi_i(r)|^2 to yield |psi_i(r')/psi_j(r')|^2
void CSUpdateBase::updateRatioMatrix(const std::vector<RealType>& ratio_i, Matrix<RealType> & Ratioij)
{
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    Ratioij[ipsi][ipsi]=1.0;

  for(int ipsi=0; ipsi< nPsi-1; ipsi++)
  {
    for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        Ratioij[ipsi][jpsi]=Ratioij[ipsi][jpsi]*ratio_i[jpsi]/ratio_i[ipsi];
        Ratioij[jpsi][ipsi]=1.0/Ratioij[ipsi][jpsi];
      }
  }
}


void CSUpdateBase::initCSWalkers(WalkerIter_t it, WalkerIter_t it_end,
                                 bool resetNorms)
{
  nPsi=Psi1.size();
  useDrift=(useDriftOption=="yes");
  if(nPsi ==0)
  {
    app_error() << "  CSUpdateBase::initCSWalkers fails. Empyty Psi/H pairs" << std::endl;
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

  UpdatePbyP=true;
  nPsi=Psi1.size();
  int nw = it_end-it;
  int iw(0);
  resizeWorkSpace(nw,W.getTotalNum());
  //ratio is for temporary holding of psi_i(...r'...)/psi_i(...r...).
  ratio.resize(nPsi,1.0);

  if(resetNorms)
    logNorm.resize(nPsi,0.0);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    avgNorm[ipsi]=std::exp(logNorm[ipsi]);
  
  for (; it != it_end; ++it)
  {
    Walker_t& awalker(**it);
    W.R=awalker.R;
    W.update(true);
    //W.loadWalker(awalker,UpdatePbyP);
//    if (awalker.MDataSet.size()!=nPsi)
//      awalker.MultDataSet.resize(nPsi);

    awalker.DataSet.clear();
    awalker.DataSet.rewind();
    awalker.registerData();
    for( int ipsi=0; ipsi<nPsi; ipsi++)
    {
      Psi1[ipsi]->registerData(W,awalker.DataSet);
    }
    awalker.DataSet.allocate();
    for( int ipsi=0; ipsi<nPsi; ipsi++)
    {
      Psi1[ipsi]->copyFromBuffer(W,awalker.DataSet);
      Psi1[ipsi]->evaluateLog(W);
      logpsi[ipsi] = Psi1[ipsi]->updateBuffer(W,awalker.DataSet,false);
      Psi1[ipsi]->G = W.G;
      *G1[ipsi] = W.G;
      Psi1[ipsi]->L = W.L;
      *L1[ipsi] = W.L;

      awalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
      awalker.Properties(ipsi,LOCALENERGY)=H1[ipsi]->evaluate(W);
      H1[ipsi]->saveProperty(awalker.getPropertyBase(ipsi));
      sumratio[ipsi]=1.0;
    }
    awalker.G=W.G;
    awalker.L=W.L;
 //   randomize(awalker);
 
    int indexij(0);
    RealType *rPtr=ratioIJ[iw++];
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
    awalker.Multiplicity=sumratio[0];
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      awalker.Properties(ipsi,UMBRELLAWEIGHT)
      = invsumratio[ipsi] =1.0/sumratio[ipsi];
      cumNorm[ipsi]+=1.0/sumratio[ipsi];
    }
  }
  #pragma omp master
  print_mem("Memory Usage after the buffer registration", app_log());
}

void CSUpdateBase::updateCSWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  int iw=0;
  while(it != it_end)
  {
    Walker_t& thisWalker(**it);
    Walker_t::WFBuffer_t& w_buffer((*it)->DataSet);
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


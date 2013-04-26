//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
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

#include "QMCWaveFunctions/Fermion/DiracDeterminantIterative.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
//extern "C"{
//#include "ILUGMRESInterface.h"
//}

namespace qmcplusplus
{

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
DiracDeterminantIterative::DiracDeterminantIterative(SPOSetBasePtr const &spos, int first):
  DiracDeterminantBase(spos,first)
{
}

///default destructor
DiracDeterminantIterative::~DiracDeterminantIterative() {}

DiracDeterminantIterative& DiracDeterminantIterative::operator=(const DiracDeterminantIterative& s)
{
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  return *this;
}


DiracDeterminantIterative::DiracDeterminantIterative(const DiracDeterminantIterative& s):
  DiracDeterminantBase(s)
{
  this->resize(s.NumPtcls,s.NumOrbitals);
}


void DiracDeterminantIterative::set(int first, int nel)
{
  APP_ABORT("DiracDeterminantIterative::set");
}

void
DiracDeterminantIterative::set_iterative(int first, int nel,double &temp_cutoff)
{
  cutoff=temp_cutoff;
  FirstIndex = first;
  resize(nel,nel);
}


void DiracDeterminantIterative::resize(int nel, int morb)
{
  particleLists.resize(nel);
  DiracDeterminantBase::resize(nel,morb);
}


void DiracDeterminantIterative::SparseToCSR(vector<int> &Arp, vector<int> &Ari,vector<double> &Arx)
{
  int systemSize=LastIndex-FirstIndex;
  int nnz_index=0;
  Arp.push_back(0);
  //    cerr<<"Particle list size is "<<particleLists.size()<<endl;
  for (int ptcl=0; ptcl<particleLists.size(); ptcl++)
  {
    for (list<pair<int,double> >::iterator orb=particleLists[ptcl].begin(); orb!=particleLists[ptcl].end(); orb++)
    {
      pair<int,double> myPair=*orb;
      int orbitalIndex=myPair.first;
      double value=myPair.second;
      //	if (abs(myValue)>=cutoff){
      Ari.push_back(orbitalIndex);
      Arx.push_back(value);
      nnz_index++;
    }
    //      }
    Arp.push_back(nnz_index);
  }
  //    cerr<<"Ari size is "<<Ari.size()<<endl;
}



DiracDeterminantBase::ValueType DiracDeterminantIterative::ratio(ParticleSet& P, int iat)
{
  //    cerr<<"Using local ratio "<<TestMe()<<endl;
  UpdateMode=ORB_PBYP_RATIO;
  WorkingIndex = iat-FirstIndex;
  assert(FirstIndex==0);
  Phi->evaluate(P, iat, psiV);
  //    cerr<<"Preparing stuff"<<endl;
  vector<int> Arp;
  vector<int> Ari;
  vector<double> Arx;
  SparseToCSR(Arp,Ari,Arx);
  vector<int> Arp2(Arp.size());
  vector<int> Ari2(Ari.size());
  vector<double> Arx2(Arx.size());
  Arp2=Arp;
  Ari2=Ari;
  Arx2=Arx;
  int particleMoved=iat;
  int systemSize=LastIndex-FirstIndex;
  int nnzUpdatedPassed=Ari.size();
  vector<double> uPassed(psiV.size());
  double detRatio_ILU=0;
  for (int i=0; i<uPassed.size(); i++)
    uPassed[i]=psiV[i];
  //    cerr<<"Calling stuff"<<systemSize<<" "<<Arp.size()<<endl;
  assert(systemSize+1==Arp.size());
  assert(Ari.size()<=nnzUpdatedPassed);
  //    cerr<<"Entering"<<endl;
  //HACK TO GET TO COMPILE    calcDeterminantILUGMRES(&particleMoved, &systemSize, &nnzUpdatedPassed, uPassed.data(), Arp.data(), Ari.data(), Arx.data(), Arp2.data(), Ari2.data(), Arx2.data(), &detRatio_ILU);
  //    int *Arp_ptr; int *Ari_ptr; double *Arx_ptr;
  //    DenseToCSR(psiM_actual,Arp_ptr,Ari_ptr,Arx_ptr);
  //    cerr<<"The size of my particle list is "<<particleLists[iat].size()<<" "<<cutoff<<endl;
  oldPtcl.clear();
  assert(iat<particleLists.size());
  particleLists[iat].swap(oldPtcl);
  for (int i=0; i<psiV.size(); i++)
  {
    if (abs(psiV(i))>=cutoff)
    {
      pair <int,double> temp(i,psiV(i));
      particleLists[iat].push_back(temp);
    }
  }
#ifdef DIRAC_USE_BLAS
  curRatio = BLAS::dot(NumOrbitals,psiM[iat-FirstIndex],&psiV[0]);
  //    cerr<<"RATIOS: "<<curRatio<<" "<<detRatio_ILU<<endl;
  return curRatio;
#else
  curRatio = DetRatio(psiM, psiV.begin(),iat-FirstIndex);
  //    cerr<<"RATIOS: "<<curRatio<<" "<<detRatio_ILU<<endl;
  return curRatio;
#endif
}



DiracDeterminantBase::ValueType DiracDeterminantIterative::ratio(ParticleSet& P, int iat,
    ParticleSet::ParticleGradient_t& dG,
    ParticleSet::ParticleLaplacian_t& dL)
{
  //    cerr<<"doing large update"<<endl;
  UpdateMode=ORB_PBYP_ALL;
  Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  WorkingIndex = iat-FirstIndex;
#ifdef DIRAC_USE_BLAS
  curRatio = BLAS::dot(NumOrbitals,psiM_temp[WorkingIndex],&psiV[0]);
#else
  curRatio= DetRatio(psiM_temp, psiV.begin(),WorkingIndex);
#endif
  if(abs(curRatio)<numeric_limits<RealType>::epsilon())
  {
    UpdateMode=ORB_PBYP_RATIO; //singularity! do not update inverse
    return 0.0;
  }
  //update psiM_temp with the row substituted
  DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,curRatio);
  //update dpsiM_temp and d2psiM_temp
  for(int j=0; j<NumOrbitals; j++)
  {
    dpsiM_temp(WorkingIndex,j)=dpsiV[j];
    d2psiM_temp(WorkingIndex,j)=d2psiV[j];
  }
  int kat=FirstIndex;
  const ValueType* restrict yptr=psiM_temp.data();
  const ValueType* restrict d2yptr=d2psiM_temp.data();
  const GradType* restrict dyptr=dpsiM_temp.data();
  for(int i=0; i<NumPtcls; i++,kat++)
  {
    //This mimics gemm with loop optimization
    GradType rv;
    ValueType lap=0.0;
    for(int j=0; j<NumOrbitals; j++,yptr++)
    {
      rv += *yptr * *dyptr++;
      lap += *yptr * *d2yptr++;
    }
    //using inline dot functions
    //GradType rv=dot(psiM_temp[i],dpsiM_temp[i],NumOrbitals);
    //ValueType lap=dot(psiM_temp[i],d2psiM_temp[i],NumOrbitals);
    //Old index: This is not pretty
    //GradType rv =psiM_temp(i,0)*dpsiM_temp(i,0);
    //ValueType lap=psiM_temp(i,0)*d2psiM_temp(i,0);
    //for(int j=1; j<NumOrbitals; j++) {
    //  rv += psiM_temp(i,j)*dpsiM_temp(i,j);
    //  lap += psiM_temp(i,j)*d2psiM_temp(i,j);
    //}
    lap -= dot(rv,rv);
    dG[kat] += rv - myG[kat];
    myG_temp[kat]=rv;
    dL[kat] += lap -myL[kat];
    myL_temp[kat]=lap;
  }
  return curRatio;
}


DiracDeterminantBase::RealType
DiracDeterminantIterative::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  return DiracDeterminantBase::evaluateLog(P,buf);
}



DiracDeterminantBase::RealType
DiracDeterminantIterative::evaluateLog(ParticleSet& P,
                                       ParticleSet::ParticleGradient_t& G,
                                       ParticleSet::ParticleLaplacian_t& L)
{
  Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
  ///I think at this point psiM has particles as the first index
  for (int ptcl=0; ptcl<psiM.extent(1); ptcl++)
  {
    particleLists[ptcl].clear();
    for (int orbital=0; orbital<psiM.extent(0); orbital++)
    {
      if (abs(psiM(orbital,ptcl))>=cutoff)
      {
        pair<int,double> temp(orbital,psiM(orbital,ptcl));
        particleLists[ptcl].push_back(temp);
      }
    }
  }
  if(NumPtcls==1)
  {
    //CurrentDet=psiM(0,0);
    ValueType det=psiM(0,0);
    ValueType y=1.0/det;
    psiM(0,0)=y;
    GradType rv = y*dpsiM(0,0);
    G(FirstIndex) += rv;
    L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
    LogValue = evaluateLogAndPhase(det,PhaseValue);
  }
  else
  {
    LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
    const ValueType* restrict yptr=psiM.data();
    const ValueType* restrict d2yptr=d2psiM.data();
    const GradType* restrict dyptr=dpsiM.data();
    for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      GradType rv;
      ValueType lap=0.0;
      for(int j=0; j<NumOrbitals; j++,yptr++)
      {
        rv += *yptr * *dyptr++;
        lap += *yptr * *d2yptr++;
      }
      G(iat) += rv;
      L(iat) += lap - dot(rv,rv);
    }
  }
  return LogValue;
}



}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3265 $   $Date: 2008-10-15 09:20:33 -0500 (Wed, 15 Oct 2008) $
 * $Id: DiracDeterminantIterative.cpp 3265 2008-10-15 14:20:33Z jnkim $
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    

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


void DiracDeterminantIterative::SparseToCSR(std::vector<int> &Arp, std::vector<int> &Ari,std::vector<double> &Arx)
{
  int systemSize=LastIndex-FirstIndex;
  int nnz_index=0;
  Arp.push_back(0);
  //    std::cerr <<"Particle list size is "<<particleLists.size()<< std::endl;
  for (int ptcl=0; ptcl<particleLists.size(); ptcl++)
  {
    for (list<std::pair<int,double> >::iterator orb=particleLists[ptcl].begin(); orb!=particleLists[ptcl].end(); orb++)
    {
      std::pair<int,double> myPair=*orb;
      int orbitalIndex=myPair.first;
      double value=myPair.second;
      //	if (std::abs(myValue)>=cutoff){
      Ari.push_back(orbitalIndex);
      Arx.push_back(value);
      nnz_index++;
    }
    //      }
    Arp.push_back(nnz_index);
  }
  //    std::cerr <<"Ari size is "<<Ari.size()<< std::endl;
}



DiracDeterminantBase::ValueType DiracDeterminantIterative::ratio(ParticleSet& P, int iat)
{
  //    std::cerr <<"Using local ratio "<<TestMe()<< std::endl;
  UpdateMode=ORB_PBYP_RATIO;
  WorkingIndex = iat-FirstIndex;
  assert(FirstIndex==0);
  Phi->evaluate(P, iat, psiV);
  //    std::cerr <<"Preparing stuff"<< std::endl;
  std::vector<int> Arp;
  std::vector<int> Ari;
  std::vector<double> Arx;
  SparseToCSR(Arp,Ari,Arx);
  std::vector<int> Arp2(Arp.size());
  std::vector<int> Ari2(Ari.size());
  std::vector<double> Arx2(Arx.size());
  Arp2=Arp;
  Ari2=Ari;
  Arx2=Arx;
  int particleMoved=iat;
  int systemSize=LastIndex-FirstIndex;
  int nnzUpdatedPassed=Ari.size();
  std::vector<double> uPassed(psiV.size());
  double detRatio_ILU=0;
  for (int i=0; i<uPassed.size(); i++)
    uPassed[i]=psiV[i];
  //    std::cerr <<"Calling stuff"<<systemSize<<" "<<Arp.size()<< std::endl;
  assert(systemSize+1==Arp.size());
  assert(Ari.size()<=nnzUpdatedPassed);
  //    std::cerr <<"Entering"<< std::endl;
  //HACK TO GET TO COMPILE    calcDeterminantILUGMRES(&particleMoved, &systemSize, &nnzUpdatedPassed, uPassed.data(), Arp.data(), Ari.data(), Arx.data(), Arp2.data(), Ari2.data(), Arx2.data(), &detRatio_ILU);
  //    int *Arp_ptr; int *Ari_ptr; double *Arx_ptr;
  //    DenseToCSR(psiM_actual,Arp_ptr,Ari_ptr,Arx_ptr);
  //    std::cerr <<"The size of my particle list is "<<particleLists[iat].size()<<" "<<cutoff<< std::endl;
  oldPtcl.clear();
  assert(iat<particleLists.size());
  particleLists[iat].swap(oldPtcl);
  for (int i=0; i<psiV.size(); i++)
  {
    if (std::abs(psiV(i))>=cutoff)
    {
      pair <int,double> temp(i,psiV(i));
      particleLists[iat].push_back(temp);
    }
  }
#ifdef DIRAC_USE_BLAS
  curRatio = BLAS::dot(NumOrbitals,psiM[iat-FirstIndex],&psiV[0]);
  //    std::cerr <<"RATIOS: "<<curRatio<<" "<<detRatio_ILU<< std::endl;
  return curRatio;
#else
  curRatio = DetRatio(psiM, psiV.begin(),iat-FirstIndex);
  //    std::cerr <<"RATIOS: "<<curRatio<<" "<<detRatio_ILU<< std::endl;
  return curRatio;
#endif
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
      if (std::abs(psiM(orbital,ptcl))>=cutoff)
      {
        std::pair<int,double> temp(orbital,psiM(orbital,ptcl));
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

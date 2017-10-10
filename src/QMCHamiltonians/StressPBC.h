//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_STRESSPBC_HAMILTONIAN_H
#define QMCPLUSPLUS_STRESSPBC_HAMILTONIAN_H
#include "coulomb_types.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "OhmmsPETE/SymTensor.h"

namespace qmcplusplus
{

struct StressPBC: public QMCHamiltonianBase, public ForceBase
{

  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
  typedef TinyVector<RealType, OHMMS_DIM> PosType;

  RealType Rcut; // parameter: radial distance within which estimator is used
  int m_exp; // parameter: exponent in polynomial fit
  int N_basis; // parameter: size of polynomial basis set
  Matrix<RealType> Sinv; // terms in fitting polynomial
  Vector<RealType> h; // terms in fitting polynomial
  Vector<RealType> c; // polynomial coefficients
  // container for short-range force estimator
  
  SymTensor<RealType, OHMMS_DIM> targetconsts;
  SymTensor<RealType, OHMMS_DIM> stress_ee_const;
  SymTensor<RealType, OHMMS_DIM> stress_eI_const;
  
  
  TrialWaveFunction& Psi;
  
  bool kcdifferent;
  RealType minkc;
  
  ///source particle set
  ParticleSet& PtclTarg;
  ParticleSet& PtclA;
  ///long-range Handler
//  LRHandlerType* AB;
  LRHandlerType* AA;
  ///locator of the distance table
  int myTableIndex;
  ///number of species of A particle set
  int NumSpeciesA;
  ///number of species of B particle set
  int NumSpeciesB;
  ///number of particles of A
  int NptclA;
  ///number of particles of B
  int NptclB;
  
  ///cutoff radius of the short-range part
  RealType myRcut;
  ///radial grid
  GridType* myGrid;
  ///Always mave a radial functor for the bare coulomb
  RadFunctorType* V0;

  int MaxGridPoints;
  
  ///number of particles per species of A
  std::vector<int> NofSpeciesA;
  ///number of particles per species of B
  std::vector<int> NofSpeciesB;
  ///Zat[iat] charge for the iat-th particle of A
  std::vector<RealType> Zat;
  ///Qat[iat] charge for the iat-th particle of B
  std::vector<RealType> Qat;
  ///Zspec[spec] charge for the spec-th species of A
  std::vector<RealType> Zspec;
  ///Qspec[spec] charge for the spec-th species of B
  std::vector<RealType> Qspec;
  ///Short-range potential for each ion
  std::vector<RadFunctorType*> Vat;
  ///Short-range potential for each species
  std::vector<RadFunctorType*> Vspec;
  
  ///On initialization, signals computation of constant stress from ion.  
  ///For cloning, forces copying of existing stress tensors.  
  bool first_time;

  ParticleSet::ParticlePos_t forces_ShortRange;
  //Constructor
  StressPBC(ParticleSet& ions, ParticleSet& elns, TrialWaveFunction& Psi, bool firsttime=true);
  //"Copy" constructor
 // StressPBC(const StressPBC& aST, ParticleSet& p, TrialWaveFunction& Psi):
//	StressPBC(aST), ForceBase(aST.PtclA, p), PtclTarg(p), Psi(Psi0)
 // {
 // }



  Return_t evaluate(ParticleSet& P);

  void InitMatrix();
  void initBreakup(ParticleSet& P);
  
  

  SymTensor<RealType,OHMMS_DIM> evaluateLR_AB(ParticleSet& P);
  SymTensor<RealType,OHMMS_DIM> evaluateSR_AB(ParticleSet& P_target);
  SymTensor<RealType,OHMMS_DIM> evaluateSR_AA(ParticleSet& P);
  SymTensor<RealType,OHMMS_DIM> evaluateLR_AA(ParticleSet& P);
  SymTensor<RealType, OHMMS_DIM> evalConsts_AB();
  SymTensor<RealType, OHMMS_DIM> evalConsts_AA(ParticleSet& P);
  
  SymTensor<RealType,OHMMS_DIM> evaluateKineticSymTensor(ParticleSet& P);
  
  Return_t g_filter(RealType r);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const
  {
    registerObservablesF(h5list,gid);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables)
  {
    addObservablesStress(plist);
  }

  void setObservables(PropertySetType& plist)
  {
    setObservablesStress(plist);
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    setParticleSetStress(plist, offset);
  }
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
//  StressPBC* makeStressClone(ParticleSet& qp, TrialWaveFunction& psi);
  bool put(xmlNodePtr cur) ;

  bool get(std::ostream& os) const
  {
    os << "Ceperley Force Estimator Hamiltonian: " << pairName;
    return true;
  }
  
  void CalculateIonIonStress()
  {
     stress_IonIon=evaluateSR_AA(PtclA)+evaluateLR_AA(PtclA)+evalConsts_AA(PtclA); //+ evaluateLR_AA(PtclA);
     stress_eI_const+=evalConsts_AB();
     stress_ee_const+=evalConsts_AA(PtclTarg);
  
     app_log()<<"\n====ion-ion stress ====\n"<<stress_IonIon<< std::endl;  
  }

};

}
#endif



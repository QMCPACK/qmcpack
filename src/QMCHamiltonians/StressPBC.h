//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_STRESSPBC_HAMILTONIAN_H
#define QMCPLUSPLUS_STRESSPBC_HAMILTONIAN_H
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


  double Rcut; // parameter: radial distance within which estimator is used
  int m_exp; // parameter: exponent in polynomial fit
  int N_basis; // parameter: size of polynomial basis set
  Matrix<RealType> Sinv; // terms in fitting polynomial
  Vector<double> h; // terms in fitting polynomial
  Vector<double> c; // polynomial coefficients
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
  vector<int> NofSpeciesA;
  ///number of particles per species of B
  vector<int> NofSpeciesB;
  ///Zat[iat] charge for the iat-th particle of A
  vector<RealType> Zat;
  ///Qat[iat] charge for the iat-th particle of B
  vector<RealType> Qat;
  ///Zspec[spec] charge for the spec-th species of A
  vector<RealType> Zspec;
  ///Qspec[spec] charge for the spec-th species of B
  vector<RealType> Qspec;
  ///Short-range potential for each ion
  vector<RadFunctorType*> Vat;
  ///Short-range potential for each species
  vector<RadFunctorType*> Vspec;
  
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

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void registerObservables(vector<observable_helper*>& h5list, hid_t gid) const
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
  
     app_log()<<"\n====ion-ion stress ====\n"<<stress_IonIon<<endl;  
  }

};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: StressPBC.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/


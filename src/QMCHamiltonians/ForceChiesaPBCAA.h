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
#ifndef QMCPLUSPLUS_FORCE_CHIESA_HAMILTONIAN_H
#define QMCPLUSPLUS_FORCE_CHIESA_HAMILTONIAN_H
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{

struct ForceChiesaPBCAA: public QMCHamiltonianBase, public ForceBase
{

  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;


  double Rcut; // parameter: radial distance within which estimator is used
  int m_exp; // parameter: exponent in polynomial fit
  int N_basis; // parameter: size of polynomial basis set
  Matrix<RealType> Sinv; // terms in fitting polynomial
  Vector<double> h; // terms in fitting polynomial
  Vector<double> c; // polynomial coefficients
  // container for short-range force estimator
  
  bool kcdifferent;
  RealType minkc;
  
  ///source particle set
  ParticleSet& PtclA;
  ///long-range Handler
  LRHandlerType* AB;
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

  bool first_time;  

  ParticleSet::ParticlePos_t forces_ShortRange;

  ForceChiesaPBCAA(ParticleSet& ions, ParticleSet& elns, bool firsttime=true);

  Return_t evaluate(ParticleSet& P);

  void InitMatrix();
  void initBreakup(ParticleSet& P);
  
  void evaluateLR(ParticleSet&);
  void evaluateSR(ParticleSet&);
  void evaluateSR_AA();
  void evaluateLR_AA();
  
  Return_t g_filter(RealType r);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void registerObservables(vector<observable_helper*>& h5list, hid_t gid) const
  {
    registerObservablesF(h5list,gid);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables);


  void setObservables(PropertySetType& plist)
  {
    QMCHamiltonianBase::setObservables(plist);
      setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    QMCHamiltonianBase::setParticlePropertyList(plist, offset);
      setParticleSetF(plist, offset);
  }


  void resetTargetParticleSet(ParticleSet& P);
  
  

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  bool put(xmlNodePtr cur) ;

  bool get(std::ostream& os) const
  {
    os << "Ceperley Force Estimator Hamiltonian: " << pairName;
    return true;
  }

};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceChiesaPBCAA.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/


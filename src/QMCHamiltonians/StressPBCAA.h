//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_STRESSPBCAA_TEMP_H
#define QMCPLUSPLUS_STRESSPBCAA_TEMP_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to StressPBCAA but uses a templated version of
 * LRHandler.
 */
struct StressPBCAA:  public QMCHamiltonianBase, public ForceBase
{

  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType       GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
  LRHandlerType* AA;
  GridType* myGrid;
  RadFunctorType* rVs;

  bool is_active;
  bool FirstTime;
  int SourceID;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NumCenters;
  SymTensor<RealType, OHMMS_DIM> myConst;
  RealType myRcut;
  string PtclRefName;
  vector<RealType> Zat,Zspec;
  vector<int> NofSpecies;
  vector<int> SpeciesID;
  
  SymTensor<RealType, OHMMS_DIM> sSR;
  SymTensor<RealType, OHMMS_DIM> sLR;

  Matrix<RealType> SR2;
  Vector<RealType> dSR;
  Vector<ComplexType> del_eikr;
  /// Flag for whether to compute forces or not
  bool ComputeForces;
//     madelung constant
  SymTensor<RealType, OHMMS_DIM> MC0;

  //single particle trace sample
 // Array<TraceReal,1>* V_sample;
 // Array<TraceReal,1>  V_const;
  ParticleSet& Ps;


  /** constructor */
  StressPBCAA(ParticleSet& ions, bool active);

  ~StressPBCAA();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  void update_source(ParticleSet& s);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

//  Return_t spevaluate(ParticleSet& P);

 // Return_t registerData(ParticleSet& P, BufferType& buffer);
//  Return_t updateBuffer(ParticleSet& P, BufferType& buffer);
 // void copyFromBuffer(ParticleSet& P, BufferType& buf);
//  void copyToBuffer(ParticleSet& P, BufferType& buf);
  Return_t evaluatePbyP(ParticleSet& P, int iat);
 // void acceptMove(int iat);

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "StressPBCAA potential: " << PtclRefName;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void initBreakup(ParticleSet& P);

  //virtual void checkout_particle_arrays(TraceManager& tm);
 //virtual void delete_particle_arrays();


  SymTensor<RealType, OHMMS_DIM> evaluateStress(){return stress;}
  SymTensor<RealType, OHMMS_DIM> evalSR(ParticleSet& P);
  SymTensor<RealType, OHMMS_DIM> evalLR(ParticleSet& P);
 // Return_t evalSRwithForces(ParticleSet& P);
//  Return_t evalLRwithForces(ParticleSet& P);
  SymTensor<RealType, OHMMS_DIM> evalConsts(bool report=true);
 // Return_t evaluateForPbyP(ParticleSet& P);

  void addObservables(PropertySetType& plist, BufferType& collectables)
  {
     addObservablesStress(plist);
  }

  void setObservables(PropertySetType& plist)
  {

      setObservablesStress(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {

     setParticleSetStress(plist, offset);
  }

};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 5976 $   $Date: 2013-09-13 13:39:44 -0500 (Fri, 13 Sep 2013) $
 * $Id: StressPBCAA.h 5976 2013-09-13 18:39:44Z jtkrogel $
 ***************************************************************************/


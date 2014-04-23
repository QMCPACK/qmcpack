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
#ifndef QMCPLUSPLUS_STRESSKINETIC_TEMP_H
#define QMCPLUSPLUS_STRESSKINETIC_TEMP_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to StressKinetic but uses a templated version of
 * LRHandler.
 */
struct StressKinetic: public QMCHamiltonianBase, public ForceBase
{

  bool FirstTime;
  int SourceID;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NumCenters;
  RealType myRcut;
  string PtclRefName;
  vector<RealType> Zat,Zspec;
  vector<int> NofSpecies;
  vector<int> SpeciesID;
  
  TrialWaveFunction& Psi;

  ParticleSet& Ps;


  /** constructor */
  StressKinetic(ParticleSet& els, TrialWaveFunction& Psi0);
  

  ~StressKinetic();
  
  Return_t evaluate(ParticleSet& P);
  Return_t evaluatePbyP(ParticleSet& P, int iat);
  
  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }
  
  SymTensor<RealType, OHMMS_DIM> evaluateKineticSymTensor(ParticleSet& P);


  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "StressKinetic potential: " << PtclRefName;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  SymTensor<RealType, OHMMS_DIM> evaluateStress(){return stress;}
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


};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 5976 $   $Date: 2013-09-13 13:39:44 -0500 (Fri, 13 Sep 2013) $
 * $Id: StressKinetic.h 5976 2013-09-13 18:39:44Z jtkrogel $
 ***************************************************************************/


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
  std::string PtclRefName;
  std::vector<RealType> Zat,Zspec;
  std::vector<int> NofSpecies;
  std::vector<int> SpeciesID;
  
  TrialWaveFunction& Psi;

  ParticleSet& Ps;


  /** constructor */
  StressKinetic(ParticleSet& els, TrialWaveFunction& Psi0);
  

  ~StressKinetic();
  
  Return_t evaluate(ParticleSet& P);
  
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
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


};

}
#endif



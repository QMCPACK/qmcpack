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
  std::string PtclRefName;
  std::vector<RealType> Zat,Zspec;
  std::vector<int> NofSpecies;
  std::vector<int> SpeciesID;
  
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

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

//  Return_t spevaluate(ParticleSet& P);

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



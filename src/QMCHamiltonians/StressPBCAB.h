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
    
    
#ifndef QMCPLUSPLUS_STRESSPBCAB_TEMP_H
#define QMCPLUSPLUS_STRESSPBCAB_TEMP_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
namespace qmcplusplus
{

/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to StressPBCAB but uses a templated version of
 * LRHandler.
 */
struct StressPBCAB: public QMCHamiltonianBase, public ForceBase
{

  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
  
  bool is_active;

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
  ///const energy after breakup
  SymTensor<RealType,OHMMS_DIM> myConst;
  ///cutoff radius of the short-range part
  RealType myRcut;
  ///radial grid
  GridType* myGrid;
  ///Always mave a radial functor for the bare coulomb
  RadFunctorType* V0;
  /// Flag for whether to compute forces or not
  bool ComputeForces;
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
  /*@{
   * @brief temporary data for pbyp evaluation
   */
  ///short-range part for the moved particle
  RealType SRtmp;
  ///long-range part for the moved particle
  RealType LRtmp;
  ///short-range per particle
 // Vector<RealType> SRpart;
  ///long-range per particle
 // Vector<RealType> LRpart;
  /*@}*/

  //This is set to true if the K_c of structure-factors are different
  bool kcdifferent;
  RealType minkc;

  //particle trace samples
 // Array<TraceReal,1>* Ve_sample;
 // Array<TraceReal,1>* Vi_sample;
 // Array<TraceReal,1>  Ve_const;
//  Array<TraceReal,1>  Vi_const;
  ParticleSet& Peln;
  ParticleSet& Pion;

  StressPBCAB(ParticleSet& ions, ParticleSet& elns, bool computeForces=false);

  ///// copy constructor
  //StressPBCAB(const StressPBCAB& c);

  ~StressPBCAB();

  void resetTargetParticleSet(ParticleSet& P);

 // virtual void checkout_particle_arrays(TraceManager& tm);
 // virtual void delete_particle_arrays();

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

 // Return_t spevaluate(ParticleSet& P);

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "StressPBCAB potential source: " << PtclA.getName() ;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void initBreakup(ParticleSet& P);

 // Return_t evalConsts_orig(bool report=true);
 // Return_t evalSR_old(ParticleSet& P);
 // Return_t evalLR_old(ParticleSet& P);
 // Return_t evalConsts_old(bool report=true);

  SymTensor<RealType, OHMMS_DIM> evaluateStress(){return stress;}
  SymTensor<RealType, OHMMS_DIM> evalSR(ParticleSet& P);
  SymTensor<RealType, OHMMS_DIM> evalLR(ParticleSet& P);
 // Return_t evalSRwithForces(ParticleSet& P);
 // Return_t evalLRwithForces(ParticleSet& P);
  SymTensor<RealType, OHMMS_DIM> evalConsts(bool report=true);
 // void add(int groupID, RadFunctorType* ppot);
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



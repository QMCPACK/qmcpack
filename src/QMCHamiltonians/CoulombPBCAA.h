//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COULOMBPBCAA_TEMP_H
#define QMCPLUSPLUS_COULOMBPBCAA_TEMP_H
#include "QMCHamiltonians/OperatorBase.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to CoulombPBCAA but uses a templated version of
 * LRHandler.
 */
struct CoulombPBCAA : public OperatorBase, public ForceBase
{
  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
  typedef LRHandlerType::mRealType mRealType;

  // using shared_ptr on AA, dAA is a compromise
  // When ion-ion is_active = false, makeClone calls the copy constructor.
  // AA, dAA are shared between clones.
  // When elec-elec is_active = true, makeClone calls the constructor
  // AA, dAA are not shared between clones and behave more like unique_ptr

  // energy-optimized
  std::shared_ptr<LRHandlerType> AA;
  GridType* myGrid;
  RadFunctorType* rVs;
  // force-optimized
  std::shared_ptr<LRHandlerType> dAA;
  GridType* myGridforce;
  RadFunctorType* rVsforce;

  bool is_active;
  bool FirstTime;
  int SourceID;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NumCenters;
  Return_t myConst;
  RealType myRcut;
  std::string PtclRefName;
  std::vector<RealType> Zat, Zspec;
  std::vector<int> NofSpecies;
  std::vector<int> SpeciesID;

  Matrix<RealType> SR2;
  Vector<RealType> dSR;
  Vector<ComplexType> del_eikr;
  /// Flag for whether to compute forces or not
  bool ComputeForces;
  //     madelung constant
  RealType MC0;

#if !defined(REMOVE_TRACEMANAGER)
  //single particle trace sample
  Array<TraceReal, 1>* V_sample;
  Array<TraceReal, 1> V_const;
#endif
  ParticleSet& Ps;


  /** constructor */
  CoulombPBCAA(ParticleSet& ref, bool active, bool computeForces = false);

  ~CoulombPBCAA();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  Return_t evaluateWithIonDerivs(ParticleSet& P,
                                 ParticleSet& ions,
                                 TrialWaveFunction& psi,
                                 ParticleSet::ParticlePos_t& hf_terms,
                                 ParticleSet::ParticlePos_t& pulay_terms);
  void update_source(ParticleSet& s);

  /** Do nothing */
  bool put(xmlNodePtr cur) { return true; }

  bool get(std::ostream& os) const
  {
    os << "CoulombPBCAA potential: " << PtclRefName;
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override;

  void initBreakup(ParticleSet& P);

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contribute_particle_quantities();
  virtual void checkout_particle_quantities(TraceManager& tm);
  Return_t evaluate_sp(ParticleSet& P); //collect
  virtual void delete_particle_quantities();
#endif

  Return_t evalSR(ParticleSet& P);
  Return_t evalLR(ParticleSet& P);
  Return_t evalSRwithForces(ParticleSet& P);
  Return_t evalLRwithForces(ParticleSet& P);
  Return_t evalConsts(bool report = true);

  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist)
  {
    OperatorBase::setObservables(plist);
    if (ComputeForces)
      setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    OperatorBase::setParticlePropertyList(plist, offset);
    if (ComputeForces)
      setParticleSetF(plist, offset);
  }

private:
  // AA table ID
  const int d_aa_ID;
};

} // namespace qmcplusplus
#endif

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COULOMBPBCAA_H
#define QMCPLUSPLUS_COULOMBPBCAA_H
#include <ResourceHandle.h>
#include "QMCHamiltonians/OperatorBase.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

template<class T>
class OneDimCubicSplineLinearGrid;

/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to CoulombPBCAA but uses a templated version of
 * LRHandler.
 */
struct CoulombPBCAA : public OperatorBase, public ForceBase
{
  using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
  using GridType       = LRCoulombSingleton::GridType;
  using RadFunctorType = LRCoulombSingleton::RadFunctorType;
  using mRealType      = LRHandlerType::mRealType;
  using OffloadSpline  = OneDimCubicSplineLinearGrid<LRCoulombSingleton::pRealType>;

  /// energy-optimized long range handle. Should be const LRHandlerType eventually
  std::shared_ptr<LRHandlerType> AA;
  /// energy-optimized short range pair potential
  std::shared_ptr<const RadFunctorType> rVs;
  /// the same as rVs but can be used inside OpenMP offload regions
  std::shared_ptr<const OffloadSpline> rVs_offload;
  /// force-optimized long range handle
  std::shared_ptr<const LRHandlerType> dAA;
  /// force-optimized short range pair potential
  std::shared_ptr<const RadFunctorType> rVsforce;

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
  std::shared_ptr<Vector<RealType, OffloadPinnedAllocator<RealType>>> Zat_offload;

  std::vector<int> NofSpecies;
  std::vector<int> SpeciesID;

  Matrix<RealType> SR2;
  Vector<RealType> dSR;
  Vector<ComplexType> del_eikr;
  /// Flag for whether to compute forces or not
  bool ComputeForces;
  /// Flag for whether to use quasi-2D Ewald
  const bool quasi2d;

#if !defined(REMOVE_TRACEMANAGER)
  //single particle trace sample
  Array<TraceReal, 1>* V_sample;
  Array<TraceReal, 1> V_const;
#endif
  ParticleSet& Ps;


  /** constructor */
  CoulombPBCAA(ParticleSet& ref, bool active, bool computeForces, bool use_offload);

  ~CoulombPBCAA() override;

  std::string getClassName() const override { return "CoulombPBCAA"; }

  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  void mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                   const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                   const RefVectorWithLeader<ParticleSet>& p_list) const override;

  /**
   * Evaluate the contribution of this component of multiple walkers per particle reporting
   * to registered listeners from Estimators.
   */
  void mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const std::vector<ListenerVector<RealType>>& listeners,
                              const std::vector<ListenerVector<RealType>>& ion_listeners) const override;


  Return_t evaluateWithIonDerivs(ParticleSet& P,
                                 ParticleSet& ions,
                                 TrialWaveFunction& psi,
                                 ParticleSet::ParticlePos& hf_terms,
                                 ParticleSet::ParticlePos& pulay_terms) override;
  void updateSource(ParticleSet& s) override;

  /** Do nothing */
  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "CoulombPBCAA potential: " << PtclRefName;
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override;

  /** Inform objects associated with this operator of per particle listeners.
   *  i.e. turnOnPerParticleSK of particleset qp.
   */
  void informOfPerParticleListener() override;

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override;
  void checkoutParticleQuantities(TraceManager& tm) override;
  Return_t evaluate_sp(ParticleSet& P); //collect
  void deleteParticleQuantities() override;
#endif

  Return_t evalSR(ParticleSet& P);

  static std::vector<Return_t> mw_evalSR_offload(const RefVectorWithLeader<OperatorBase>& o_list,
                                                 const RefVectorWithLeader<ParticleSet>& p_list);

  Return_t evalLR(ParticleSet& P);
  Return_t evalSRwithForces(ParticleSet& P);
  Return_t evalLRwithForces(ParticleSet& P);
  Return_t evalConsts(bool report = true);

  void addObservables(PropertySetType& plist, BufferType& collectables) override;

  void setObservables(PropertySetType& plist) override
  {
    OperatorBase::setObservables(plist);
    if (ComputeForces)
      setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset) override
  {
    OperatorBase::setParticlePropertyList(plist, offset);
    if (ComputeForces)
      setParticleSetF(plist, offset);
  }

  /** initialize a shared resource and hand it to a collection
   */
  void createResource(ResourceCollection& collection) const override;

  /** acquire a shared resource from a collection
   */
  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

  /** return a shared resource to a collection
   */
  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

  RealType get_madelung_constant() const { return madelung_constant_; }

private:
  RealType madelung_constant_;

  /// if true use offload
  const bool use_offload_;
  /// AA table ID
  const int d_aa_ID;
  /// Timer for long range
  NewTimer& evalLR_timer_;
  /// Timer for long range
  NewTimer& evalSR_timer_;
  /// Timer for offload part
  NewTimer& offload_timer_;

  /// multiwalker shared resource
  struct CoulombPBCAAMultiWalkerResource;
  ResourceHandle<CoulombPBCAAMultiWalkerResource> mw_res_handle_;

  /** constructor code factored out
   */
  void initBreakup(ParticleSet& P);

  /** Compute the const part of the per particle coulomb self interaction potential.
   *  \param[out]  pp_consts   constant values for the particles self interaction
   */
  void evalPerParticleConsts(Vector<RealType>& pp_consts) const;
};

} // namespace qmcplusplus
#endif

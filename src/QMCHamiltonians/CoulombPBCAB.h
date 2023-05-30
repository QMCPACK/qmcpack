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


#ifndef QMCPLUSPLUS_COULOMBPBCAB_H
#define QMCPLUSPLUS_COULOMBPBCAB_H
#include <ResourceHandle.h>
#include "QMCHamiltonians/OperatorBase.h"
#include "QMCHamiltonians/ForceBase.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to CoulombPBCAB but uses a templated version of
 * LRHandler.
 */
class CoulombPBCAB : public OperatorBase, public ForceBase
{
public:
  using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
  using GridType       = LRCoulombSingleton::GridType;
  using RadFunctorType = LRCoulombSingleton::RadFunctorType;
  using mRealType      = LRHandlerType::mRealType;

  ///long-range Handler. Should be const LRHandlerType eventually
  std::shared_ptr<LRHandlerType> AB;
  ///long-range derivative handler
  std::shared_ptr<const LRHandlerType> dAB;
  ///locator of the distance table
  const int myTableIndex;
  ///number of species of A particle set
  int NumSpeciesA;
  ///number of species of B particle set
  int NumSpeciesB;
  ///number of particles of A
  int NptclA;
  ///number of particles of B
  int NptclB;
  ///const energy after breakup
  Return_t myConst;
  ///cutoff radius of the short-range part
  RealType myRcut;
  ///radial grid
  std::shared_ptr<GridType> myGrid;
  ///Always mave a radial functor for the bare coulomb
  std::shared_ptr<const RadFunctorType> V0;
  ///Radial functor for bare coulomb, optimized for forces
  std::shared_ptr<const RadFunctorType> fV0;
  ///Radial functor for derivative of bare coulomb, optimized for forces
  std::shared_ptr<const RadFunctorType> dfV0;
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
  std::vector<const RadFunctorType*> Vat;
  ///Short-range potential for each species
  std::vector<std::shared_ptr<RadFunctorType>> Vspec;
  ///Short-range potential (r*V) and potential derivative d/dr(rV) derivative for each ion
  ///Required for force evaluations.
  std::vector<const RadFunctorType*> fVat;
  std::vector<const RadFunctorType*> fdVat;
  ////Short-range potential (r*V) and potential derivative d/dr(rV) derivative for each species
  std::vector<std::shared_ptr<const RadFunctorType>> fVspec;
  std::vector<std::shared_ptr<const RadFunctorType>> fdVspec;
  /*@{
   * @brief temporary data for pbyp evaluation
   */
  ///short-range part for the moved particle
  RealType SRtmp;
  ///long-range part for the moved particle
  RealType LRtmp;
  ///short-range per particle
  Vector<RealType> SRpart;
  ///long-range per particle
  Vector<RealType> LRpart;
  /*@}*/

  //This is set to true if the K_c of structure-factors are different
  bool kcdifferent;
  RealType minkc;

#if !defined(REMOVE_TRACEMANAGER)
  //particle trace samples
  Array<TraceReal, 1>* Ve_sample;
  Array<TraceReal, 1>* Vi_sample;
  Array<TraceReal, 1> Ve_const;
  Array<TraceReal, 1> Vi_const;
#endif
  // \todo Coulomb class is walker agnositic, it should not record a particular electron particle set.
  // kept for the trace manager. Delete this particle set reference when support for TraceManager is permanently
  // removed which should coincide with the removal of the legacy drivers.
  ParticleSet& Peln;

  CoulombPBCAB(ParticleSet& ions, ParticleSet& elns, bool computeForces = false);

  ///// copy constructor
  //CoulombPBCAB(const CoulombPBCAB& c);

  ~CoulombPBCAB() override;

  void resetTargetParticleSet(ParticleSet& P) override;

  std::string getClassName() const override { return "CoulombPBCAB"; }

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override;
  void checkoutParticleQuantities(TraceManager& tm) override;
  Return_t evaluate_sp(ParticleSet& P); //collect
  void deleteParticleQuantities() override;
#endif


  Return_t evaluate(ParticleSet& P) override;

  /** Evaluate the contribution of this component of multiple walkers per particle reporting
   *  to registered listeners from Estimators.
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

  /** Do nothing */
  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "CoulombPBCAB potential source: " << pset_ions_.getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override;

  ///Computes the short-range contribution to the coulomb energy.
  Return_t evalSR(ParticleSet& P);
  ///Computes the long-range contribution to the coulomb energy.
  Return_t evalLR(ParticleSet& P);
  ///Computes the short-range contribution to the coulomb energy and forces.
  Return_t evalSRwithForces(ParticleSet& P);
  ///Computes the long-range contribution to the coulomb energy and forces.
  Return_t evalLRwithForces(ParticleSet& P);
  ///Evaluates madelung and background contributions to total energy.
  Return_t evalConsts(const ParticleSet& P, bool report = true);
  ///Adds a local pseudopotential channel "ppot" to all source species of type "groupID".
  void add(int groupID, std::unique_ptr<RadFunctorType>&& ppot);

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

  const ParticleSet& getSourcePSet() const { return pset_ions_; }

  /** initialize a shared resource and hand it to a collection
   */
  void createResource(ResourceCollection& collection) const override;

  /** acquire a shared resource from a collection
   */
  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

  /** return a shared resource to a collection
   */
  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

  /** Call to inform objects associated with this operator of per particle listeners.
   *  should be called before createResources
   */
  void informOfPerParticleListener() override;

protected:
  /** Creates the long-range handlers, then splines and stores it by particle and species for quick evaluation.
   *  this is just constructor code factored out.
   *  It is called by the derived class CoulombPBCAB_CUDA
   */
  void initBreakup(ParticleSet& P);

private:
  ///source particle set
  ParticleSet& pset_ions_;

  struct CoulombPBCABMultiWalkerResource;
  ResourceHandle<CoulombPBCABMultiWalkerResource> mw_res_handle_;

  /** Compute the const part of the per particle coulomb AB potential.
   *  \param[out]  pp_consts_src   constant values for the source particles aka ions aka A   
   *  \param[out]  pp_consts_trg   constant values for the target particles aka electrons aka B
   */
  void evalPerParticleConsts(Vector<RealType>& pp_consts_src, Vector<RealType>& pp_consts_trg) const;
};

} // namespace qmcplusplus
#endif

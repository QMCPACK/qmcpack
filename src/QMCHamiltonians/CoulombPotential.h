//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COULOMBPOTENTIAL_H
#define QMCPLUSPLUS_COULOMBPOTENTIAL_H
#include <ParticleSet.h>
#include <DistanceTable.h>
#include <MCWalkerConfiguration.h>
#include "ForceBase.h"
#include "OperatorBase.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

/** CoulombPotential
 * @tparam T type of the elementary data
 *
 * Hamiltonian operator for the Coulomb interaction for both AA and AB type for open systems.
 */

class CoulombPotential : public OperatorDependsOnlyOnParticleSet, public ForceBase
{
public:
  struct CoulombPotentialMultiWalkerResource;

private:
  ///source particle set
  ParticleSet& Pa;
  ///target particle set
  ParticleSet& Pb;
  ///distance table index
  const int myTableIndex;
  ///true if the table is AA
  const bool is_AA;
  ///true, if CoulombAA for quantum particleset
  bool is_active;
  ///number of centers
  int nCenters;
#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace samples
  Array<TraceReal, 1>* Va_sample;
  Array<TraceReal, 1>* Vb_sample;
#endif

  /// Flag for whether to compute forces or not
  bool ComputeForces;

public:
  /** constructor for AA
   * @param s source particleset
   * @param active if true, new Value is computed whenver evaluate is used.
   * @param computeForces if true, computes forces between inactive species
   *
   * in practice active is true for electrons and false for ions.
   * we have no choice but to use this later to determine which listeners AA Coulomb potentials should
   * publish to.
   */
  CoulombPotential(ParticleSet& s, bool active, bool computeForces, bool copy = false);

  /** constructor for AB
   * @param s source particleset
   * @param t target particleset
   * @param active if true, new Value is computed whenver evaluate is used.
   * @param ComputeForces is not implemented for AB
   */
  CoulombPotential(ParticleSet& s, ParticleSet& t, bool active, bool copy = false);

  std::string getClassName() const override;

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override;
  void checkoutParticleQuantities(TraceManager& tm) override;
  void deleteParticleQuantities() override;
#endif

  void addObservables(PropertySetType& plist, BufferType& collectables) override;

  /** evaluate AA-type interactions */
  Return_t evaluateAA(const DistanceTableAA& d, const ParticleScalar* restrict Z);

  /** evaluate AA-type forces */
  void evaluateAAForces(const DistanceTableAA& d, const ParticleScalar* restrict Z);

  /** JNKIM: Need to check the precision */
  Return_t evaluateAB(const DistanceTableAB& d, const ParticleScalar* restrict Za, const ParticleScalar* restrict Zb);

  static Return_t evaluate_spAA(const DistanceTableAA& d,
                                const ParticleScalar* restrict Z,
                                Vector<RealType>& ve_samples,
                                const std::vector<ListenerVector<RealType>>& listeners);
  static Return_t evaluate_spAB(const DistanceTableAB& d,
                                const ParticleScalar* restrict Za,
                                const ParticleScalar* restrict Zb,
                                Vector<RealType>& ve_samples,
                                Vector<RealType>& vi_samples,
                                const std::vector<ListenerVector<RealType>>& listeners,
                                const std::vector<ListenerVector<RealType>>& ion_listeners);

#if !defined(REMOVE_TRACEMANAGER)
  /** evaluate AA-type interactions */
  Return_t evaluate_spAA(const DistanceTableAA& d, const ParticleScalar* restrict Z);
  Return_t evaluate_spAB(const DistanceTableAB& d,
                         const ParticleScalar* restrict Za,
                         const ParticleScalar* restrict Zb);
#endif
  ~CoulombPotential() override {};

  void updateSource(ParticleSet& s) override;
  Return_t evaluate(ParticleSet& P) override;

  void mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const std::vector<ListenerVector<RealType>>& listeners,
                              const std::vector<ListenerVector<RealType>>& ion_listeners) const override;

  void evaluateIonDerivs(ParticleSet& P,
                         ParticleSet& ions,
                         TrialWaveFunction& psi,
                         ParticleSet::ParticlePos& hf_terms,
                         ParticleSet::ParticlePos& pulay_terms) override;

  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override;

  void setObservables(PropertySetType& plist) override;

  void setParticlePropertyList(PropertySetType& plist, int offset) override;

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp) override;

private:
  ResourceHandle<CoulombPotentialMultiWalkerResource> mw_res_handle_;
};

} // namespace qmcplusplus
#endif

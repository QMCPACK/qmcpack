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


#ifndef QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#define QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_H
#include "Configuration.h"
#include "QMCHamiltonians/ForceBase.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "NeighborListsForPseudo.h"
#include "type_traits/OptionalRef.hpp"

namespace qmcplusplus
{
class NonLocalECPComponent;
template<typename T>
struct NLPPJob;

namespace testing
{
class TestNonLocalECPotential;
}

/** @ingroup hamiltonian
 * \brief Evaluate the semi local potentials
 */
class NonLocalECPotential : public OperatorBase, public ForceBase
{
  using Real = QMCTraits::RealType;

  struct NonLocalECPotentialMultiWalkerResource;

public:
  NonLocalECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi, bool enable_DLA, bool use_VP);
  NonLocalECPotential(const NonLocalECPotential& nlpp, ParticleSet& els, TrialWaveFunction& psi);
  ~NonLocalECPotential() override;

  bool dependsOnWaveFunction() const override { return true; }
  std::string getClassName() const override { return "NonLocalECPotential"; }
  void resetTargetParticleSet(ParticleSet& P) override;

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override;
  void checkoutParticleQuantities(TraceManager& tm) override;
  void deleteParticleQuantities() override;
#endif

  Return_t evaluate(ParticleSet& P) override;
  Return_t evaluateDeterministic(ParticleSet& P) override;
  void mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                   const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                   const RefVectorWithLeader<ParticleSet>& p_list) const override;

  Return_t evaluateWithToperator(ParticleSet& P) override;

  void mw_evaluateWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                const RefVectorWithLeader<ParticleSet>& p_list) const override;

  void mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const std::vector<ListenerVector<Real>>& listeners,
                              const std::vector<ListenerVector<Real>>& listeners_ions) const override;

  void mw_evaluatePerParticleWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                           const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           const std::vector<ListenerVector<Real>>& listeners,
                                           const std::vector<ListenerVector<Real>>& listeners_ions) const override;

  void evaluateIonDerivs(ParticleSet& P,
                         ParticleSet& ions,
                         TrialWaveFunction& psi,
                         ParticleSet::ParticlePos& hf_terms,
                         ParticleSet::ParticlePos& pulay_terms) override;

  void evaluateOneBodyOpMatrix(ParticleSet& P, const TWFFastDerivWrapper& psi, std::vector<ValueMatrix>& B) override;

  void evaluateOneBodyOpMatrixForceDeriv(ParticleSet& P,
                                         ParticleSet& source,
                                         const TWFFastDerivWrapper& psi,
                                         const int iat,
                                         std::vector<std::vector<ValueMatrix>>& Bforce) override;


  /** make non local moves with particle-by-particle moves
   * @param P particle set
   * @return the number of accepted moves
   */
  int makeNonLocalMovesPbyP(ParticleSet& P, NonLocalTOperator& move_op) override;

  Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       const Vector<ValueType>& dlogpsi,
                                       Vector<ValueType>& dhpsioverpsi) override;

  /** Do nothing */
  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "NonLocalECPotential: " << IonConfig.getName();
    return true;
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

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override;

  void addComponent(int groupID, std::unique_ptr<NonLocalECPComponent>&& pp);

  /** set the internal RNG pointer as the given pointer
   * @param rng input RNG pointer
   */
  void setRandomGenerator(RandomBase<FullPrecRealType>* rng) override { myRNG = rng; }

protected:
  /** the actual implementation for batched walkers, used by mw_evaluate, mw_evaluateWithToperator
   *  mw_evaluatePerPaticleWithToperator
   * @param o_list     the list of NonLocalECPotential in a walker batch
   * @param wf_list    the list of TrialWaveFunction in a walker batch
   * @param p_list     the list of ParticleSet in a walker batch
   * @param compute_txy_all whether to compute Txy for all the electrons affected by NLPP
   * @param listeners  optional listeners which allow per particle and reduced to share impl
   */
  static void mw_evaluateImpl(const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              bool compute_txy_all,
                              std::optional<ListenerOption<Real>> listeners,
                              bool keepGrid = false);

  ///random number generator
  RandomBase<FullPrecRealType>* myRNG;
  ///the set of local-potentials (one for each ion)
  std::vector<NonLocalECPComponent*> PP;
  ///unique NonLocalECPComponent to remove
  std::vector<std::unique_ptr<NonLocalECPComponent>> PPset;
  ///reference to the center ion
  ParticleSet& IonConfig;
  ///target TrialWaveFunction
  TrialWaveFunction& Psi;
  ///true, determinant localization approximation(DLA) is enabled
  bool use_DLA;

private:
  ///virtual particle set
  const std::unique_ptr<VirtualParticleSet> vp_;
  ///number of ions
  int NumIons;
  ///index of distance table for the ion-el pair
  int myTableIndex;
  ///reference to the electrons
  ParticleSet& Peln;
  ///neighbor lists for marking electrons touched by T-moves
  NeighborListsForPseudo neighbor_lists;
  ///ture if an electron is affected by other electrons moved by T-moves
  std::vector<bool> elecTMAffected;
  ///Pulay force vector
  ParticleSet::ParticlePos PulayTerm;
  /// Tmove data collected for all the electrons
  std::vector<NonLocalData> tmove_xy_all_;
#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace samples

  Array<TraceReal, 1>* Ve_sample;
  Array<TraceReal, 1>* Vi_sample;
#endif
  ///NLPP job list of ion-electron pairs by spin group
  std::vector<std::vector<NLPPJob<Real>>> nlpp_jobs;
  /// mult walker shared resource
  ResourceHandle<NonLocalECPotentialMultiWalkerResource> mw_res_handle_;

  /** the actual implementation, used by evaluate and evaluateWithToperator
   * @param P particle set
   * @param compute_txy_all whether to compute Txy for all the electrons affected by NLPP
   * @param keepGrid.  If true, does not randomize the quadrature grid before evaluation.  
   */
  void evaluateImpl(ParticleSet& P, bool compute_txy_all, bool keepGrid = false);

  /** compute the T move transition probability for a given electron
   * member variable nonLocalOps.Txy is updated
   * @param P particle set
   * @param ref_elec reference electron id
   * @param tmove_xy off-diagonal terms for one electron.
   */
  void computeOneElectronTxy(ParticleSet& P, const int ref_elec, std::vector<NonLocalData>& tmove_xy);

  friend class testing::TestNonLocalECPotential;
};
} // namespace qmcplusplus
#endif

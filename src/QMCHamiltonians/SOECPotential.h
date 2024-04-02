//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SO_ECPOTENTIAL_H
#define QMCPLUSPLUS_SO_ECPOTENTIAL_H

#include "QMCHamiltonians/SOECPComponent.h"
#include "Particle/NeighborLists.h"

namespace qmcplusplus
{
template<typename T>
struct NLPPJob;

namespace testing
{
class TestSOECPotential;
}

class SOECPotential : public OperatorBase
{
  struct SOECPotentialMultiWalkerResource;
  using Real = QMCTraits::RealType;

public:
  SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);
  ~SOECPotential() override;

  bool dependsOnWaveFunction() const override { return true; }
  std::string getClassName() const override { return "SOECPotential"; }
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;
  Return_t evaluateDeterministic(ParticleSet& P) override;

  Return_t evaluateValueAndDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       const Vector<ValueType>& dlogpsi,
                                       Vector<ValueType>& dhpsioverpsi) override;

  void mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                   const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                   const RefVectorWithLeader<ParticleSet>& p_list) const override;

  void mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const std::vector<ListenerVector<Real>>& listeners,
                              const std::vector<ListenerVector<Real>>& listeners_ions) const override;

  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "SOECPotential: " << ion_config_.getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  void addComponent(int groupID, std::unique_ptr<SOECPComponent>&& pp);

  void setRandomGenerator(RandomBase<FullPrecRealType>* rng) override { my_rng_ = rng; }

  //initialize shared resource and hand to collection
  void createResource(ResourceCollection& collection) const override;

  //acquire a shared resource
  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

  //return a shared resource to a collection
  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<OperatorBase>& o_list) const override;

protected:
  RandomBase<FullPrecRealType>* my_rng_;
  std::vector<SOECPComponent*> pp_;
  std::vector<std::unique_ptr<SOECPComponent>> ppset_;
  ParticleSet& ion_config_;
  TrialWaveFunction& psi_;
  static void mw_evaluateImpl(const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              std::optional<ListenerOption<Real>> listeners,
                              bool keep_grid = false);

private:
  ///number of ions
  int num_ions_;
  ///index of distance table for ion-el pair
  int my_table_index_;
  ///reference to the electrons
  ParticleSet& peln_;
  ///neighborlist of electrons
  NeighborLists elec_neighbor_ions_;
  ///neighborlist of ions
  NeighborLists ion_neighbor_elecs_;
  //job list for evaluation
  std::vector<std::vector<NLPPJob<RealType>>> sopp_jobs_;
  //multi walker resource
  ResourceHandle<SOECPotentialMultiWalkerResource> mw_res_handle_;

  void evaluateImpl(ParticleSet& elec, bool keep_grid = false);

  //for testing
  friend class testing::TestSOECPotential;
};

} // namespace qmcplusplus

#endif

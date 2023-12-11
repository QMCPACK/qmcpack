//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "CostFunctionCrowdData.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{
CostFunctionCrowdData::CostFunctionCrowdData(int crowd_size,
                                             ParticleSet& P,
                                             TrialWaveFunction& Psi,
                                             QMCHamiltonian& H,
                                             RandomBase<FullPrecRealType>& Rng)
    : h0_res_("h0 resource"), e0_(0.0), e2_(0.0), wgt_(0.0), wgt2_(0.0)
{
  P.createResource(driverwalker_resource_collection_.pset_res);
  Psi.createResource(driverwalker_resource_collection_.twf_res);
  H.createResource(driverwalker_resource_collection_.ham_res);

  log_psi_fixed_.resize(crowd_size);
  log_psi_opt_.resize(crowd_size);

  wf_ptr_list_.resize(crowd_size);
  p_ptr_list_.resize(crowd_size);
  h_ptr_list_.resize(crowd_size);
  h0_ptr_list_.resize(crowd_size);

  rng_ptr_list_.resize(crowd_size);

  // build a temporary H_KE for later calling makeClone
  // need makeClone to setup internal my_index_ of a new copy.
  const auto components = H.getTWFDependentComponents();
  QMCHamiltonian H_KE("h_free");
  for (OperatorBase& component : components)
    H_KE.addOperator(component.makeClone(P, Psi), component.getName());
  H_KE.createResource(h0_res_);

  for (int ib = 0; ib < crowd_size; ib++)
  {
    p_ptr_list_[ib] = std::make_unique<ParticleSet>(P);
    auto& pCopy     = *p_ptr_list_[ib];

    wf_ptr_list_[ib] = Psi.makeClone(pCopy);
    auto& psiCopy    = *wf_ptr_list_[ib];

    h_ptr_list_[ib]  = H.makeClone(pCopy, psiCopy);
    h0_ptr_list_[ib] = H_KE.makeClone(pCopy, psiCopy);

    rng_ptr_list_[ib] = Rng.makeClone();
    h_ptr_list_[ib]->setRandomGenerator(rng_ptr_list_[ib].get());
    h0_ptr_list_[ib]->setRandomGenerator(rng_ptr_list_[ib].get());
  }
  rng_save_ptr_ = Rng.makeClone();
}

RefVector<ParticleSet> CostFunctionCrowdData::get_p_list(int len)
{
  return convertUPtrToRefVectorSubset(p_ptr_list_, 0, len);
}

RefVector<TrialWaveFunction> CostFunctionCrowdData::get_wf_list(int len)
{
  return convertUPtrToRefVectorSubset(wf_ptr_list_, 0, len);
}

RefVector<QMCHamiltonian> CostFunctionCrowdData::get_h_list(int len)
{
  return convertUPtrToRefVectorSubset(h_ptr_list_, 0, len);
}

RefVector<QMCHamiltonian> CostFunctionCrowdData::get_h0_list(int len)
{
  return convertUPtrToRefVectorSubset(h0_ptr_list_, 0, len);
}

void CostFunctionCrowdData::zero_log_psi()
{
  std::fill(log_psi_opt_.begin(), log_psi_opt_.end(), 0.0);
  std::fill(log_psi_fixed_.begin(), log_psi_fixed_.end(), 0.0);
}

} // namespace qmcplusplus

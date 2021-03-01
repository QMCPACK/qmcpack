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
                                             QMCHamiltonian& H_KE_Node,
                                             RandomGenerator_t& Rng)
    : e0_(0.0), e2_(0.0), wgt_(0.0), wgt2_(0.0)
{
  log_psi_fixed_.resize(crowd_size);
  log_psi_opt_.resize(crowd_size);

  wf_ptr_list_.resize(crowd_size);
  p_ptr_list_.resize(crowd_size);
  h_ptr_list_.resize(crowd_size);
  h0_ptr_list_.resize(crowd_size);

  rng_ptr_list_.resize(crowd_size);

  for (int ib = 0; ib < crowd_size; ib++)
  {
    ParticleSet* pCopy         = new ParticleSet(P);
    TrialWaveFunction* psiCopy = Psi.makeClone(*pCopy);

    p_ptr_list_[ib].reset(pCopy);
    wf_ptr_list_[ib].reset(psiCopy);
    h_ptr_list_[ib].reset(H.makeClone(*pCopy, *psiCopy));
    h0_ptr_list_[ib].reset(H_KE_Node.makeClone(*pCopy, *psiCopy));

    rng_ptr_list_[ib] = std::make_unique<RandomGenerator_t>(Rng);
    h_ptr_list_[ib]->setRandomGenerator(rng_ptr_list_[ib].get());
    h0_ptr_list_[ib]->setRandomGenerator(rng_ptr_list_[ib].get());

    rng_save_ptr_ = std::make_unique<RandomGenerator_t>(Rng);
  }
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

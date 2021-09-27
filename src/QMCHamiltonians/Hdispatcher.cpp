//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "Hdispatcher.h"
#include <cassert>
#include "QMCHamiltonian.h"

namespace qmcplusplus
{
Hdispatcher::Hdispatcher(bool use_batch) : use_batch_(use_batch) {}


std::vector<QMCHamiltonian::FullPrecRealType> Hdispatcher::flex_evaluate(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list) const
{
  assert(ham_list.size() == p_list.size());
  if (use_batch_)
    return QMCHamiltonian::mw_evaluate(ham_list, wf_list, p_list);
  else
  {
    std::vector<FullPrecRealType> local_energies(ham_list.size());
    for (size_t iw = 0; iw < ham_list.size(); iw++)
      local_energies[iw] = ham_list[iw].evaluate(p_list[iw]);
    return local_energies;
  }
}

std::vector<QMCHamiltonian::FullPrecRealType> Hdispatcher::flex_evaluateWithToperator(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list) const
{
  assert(ham_list.size() == p_list.size());
  if (use_batch_)
    return QMCHamiltonian::mw_evaluateWithToperator(ham_list, wf_list, p_list);
  else
  {
    std::vector<FullPrecRealType> local_energies(ham_list.size());
    for (size_t iw = 0; iw < ham_list.size(); iw++)
      local_energies[iw] = ham_list[iw].evaluateWithToperator(p_list[iw]);
    return local_energies;
  }
}

std::vector<int> Hdispatcher::flex_makeNonLocalMoves(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                                     const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                     const RefVectorWithLeader<ParticleSet>& p_list) const
{
  assert(ham_list.size() == p_list.size());
  if (use_batch_)
    return QMCHamiltonian::mw_makeNonLocalMoves(ham_list, wf_list, p_list);
  else
  {
    std::vector<int> num_accepts(ham_list.size());
    for (size_t iw = 0; iw < ham_list.size(); iw++)
      num_accepts[iw] = ham_list[iw].makeNonLocalMoves(p_list[iw]);
    return num_accepts;
  }
}

} // namespace qmcplusplus

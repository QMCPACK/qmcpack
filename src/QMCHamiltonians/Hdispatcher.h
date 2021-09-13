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


#ifndef QMCPLUSPLUS_HDISPATCH_H
#define QMCPLUSPLUS_HDISPATCH_H

#include "QMCHamiltonian.h"

namespace qmcplusplus
{
/** Wrappers for dispatching to QMCHamiltonian single walker APIs or mw_ APIs.
 * This should be only used by QMC drivers.
 * member function names must match mw_ APIs in QMCHamiltonian
 */
class Hdispatcher
{
public:
  using FullPrecRealType = QMCHamiltonian::FullPrecRealType;

  Hdispatcher(bool use_batch);

  std::vector<FullPrecRealType> flex_evaluate(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                              const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                              const RefVectorWithLeader<ParticleSet>& p_list) const;

  std::vector<FullPrecRealType> flex_evaluateWithToperator(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                                           const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                           const RefVectorWithLeader<ParticleSet>& p_list) const;

  std::vector<int> flex_makeNonLocalMoves(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                          const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list) const;

private:
  bool use_batch_;
};
} // namespace qmcplusplus

#endif

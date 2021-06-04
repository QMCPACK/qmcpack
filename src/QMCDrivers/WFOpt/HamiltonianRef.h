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


#ifndef QMCPLUSPLUS_HAMILTONIANREF_H
#define QMCPLUSPLUS_HAMILTONIANREF_H

#include <OperatorBase.h>
#ifdef QMC_CUDA
class MCWalkerConfiguration;
#endif

namespace qmcplusplus
{
/* Helper class to handle references to H components.
 */
class HamiltonianRef
{
public:
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  using Walker_t         = OperatorBase::Walker_t;
  using ValueType        = OperatorBase::ValueType;
  using RealType         = OperatorBase::RealType;

  /// record operator reference
  void addOperator(OperatorBase& op);

  /// the same evaluateValueAndDerivatives as QMCHamiltonian
  FullPrecRealType evaluateValueAndDerivatives(ParticleSet& P,
                                               const opt_variables_type& optvars,
                                               std::vector<ValueType>& dlogpsi,
                                               std::vector<ValueType>& dhpsioverpsi,
                                               bool compute_deriv);

  /// the same evaluate as QMCHamiltonian
  FullPrecRealType evaluate(ParticleSet& P);

#ifdef QMC_CUDA
  /// the same evaluate as QMCHamiltonian
  void evaluate(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy);
#endif

  int addObservables(ParticleSet& P);

  int size() const { return Hrefs_.size(); }

private:
  /// collected references
  RefVector<OperatorBase> Hrefs_;
};

} // namespace qmcplusplus

#endif

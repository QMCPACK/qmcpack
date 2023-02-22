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

  HamiltonianRef(const RefVector<OperatorBase>);

  /// the same evaluateValueAndDerivatives as QMCHamiltonian
  FullPrecRealType evaluateValueAndDerivatives(ParticleSet& P,
                                               const opt_variables_type& optvars,
                                               Vector<ValueType>& dlogpsi,
                                               Vector<ValueType>& dhpsioverpsi);

  /// the same evaluate as QMCHamiltonian
  FullPrecRealType evaluate(ParticleSet& P);

  int size() const { return Hrefs_.size(); }

private:
  /// collected references
  const RefVector<OperatorBase> Hrefs_;
};

} // namespace qmcplusplus

#endif

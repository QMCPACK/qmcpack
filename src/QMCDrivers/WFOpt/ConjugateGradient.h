//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CONGJUGATEGRADIENT_H
#define QMCPLUSPLUS_CONGJUGATEGRADIENT_H


#include <Configuration.h>
#include <NewTimer.h>

namespace qmcplusplus
{
class QMCCostFunctionBase;

class ConjugateGradient
{
  using Real = QMCTraits::RealType;

public:
  ConjugateGradient(Real thr = 1e-6, Real regularization = 0);

  //Solve linear system Ax = b, i.e. x = Ainv * b using conjugate gradient.
  //If ill-conditioned, solves x = (Ainv + lambda * I) * b where lambda is a small regularization
  int run(QMCCostFunctionBase& optTarget, const std::vector<Real>& bvec, std::vector<Real>& solution);

  Real getNonLinearRescale(const QMCCostFunctionBase& optTarget) const;

private:
  //convergence threshold
  Real threshold_;
  //for ill-conditioned matrices, add small offset to the diagonal
  Real regularization_;
  //to be used for nonlinear rescale
  std::vector<Real> Axsol_;
  std::vector<Real> xsol_;
  NewTimer& conjugate_gradient_timer_;

  // obtain the range of non-linear parameters
  void getNonLinearRange(int& first, int& last, const QMCCostFunctionBase& optTarget) const;
};

} // namespace qmcplusplus
#endif

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

private:

  Real threshold_;
  Real regularization_;
};

}
#endif

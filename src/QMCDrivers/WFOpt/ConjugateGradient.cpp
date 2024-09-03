//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "ConjugateGradient.h"
#include "QMCCostFunctionBase.h"

namespace qmcplusplus
{

ConjugateGradient::ConjugateGradient(Real threshold, Real regularization)
    : threshold_(threshold), regularization_(regularization)
{}

int ConjugateGradient::run(QMCCostFunctionBase& optTarget, const std::vector<Real>& bvec, std::vector<Real>& solution)
{
  int numParams = optTarget.getNumParams();
  int kmax = numParams;
  int k    = 0;

  std::vector<Real> rk(numParams, 0);
  std::vector<Real> rkp1(numParams, 0);
  std::vector<Real> pk(numParams, 0);
  std::vector<Real> pkp1(numParams, 0);
  std::vector<Real> xk(numParams, 0);
  std::vector<Real> xkp1(numParams, 0);
  std::vector<Real> Apk(numParams, 0);
  solution.resize(numParams);


  //Basic implementation of conjugate gradient algorithm...coming directly from wikipedia
  Real dk = 0;
  //initial guess is zero, so S*x_0 = 0
  for (int i = 0; i < numParams; i++)
  {
    rk[i] = bvec[i];
    dk += rk[i] * rk[i];
  }
  Real eps       = dk / numParams * threshold_;
  pk             = rk;
  bool converged = false;
  while (!converged)
  {
    //This function avoids building the full Sij matrix by calculating the matrix vector
    //product of the overlap matrix * the current parameter vector.
    //This is the "accelerated stochastic reconfiguration" from https://doi.org/10.1103/PhysRevB.85.045103
    //Basically doing Eqn. (7)
    optTarget.calcOvlParmVec(pk, Apk);
    for (int pm = 0; pm < numParams; pm++)
      Apk[pm] += regularization_ * pk[pm];
    Real denom = 0;
    for (int i = 0; i < numParams; i++)
      denom += pk[i] * Apk[i];
    Real ak = dk / denom;

    Real dkp1 = 0;
    for (int i = 0; i < numParams; i++)
    {
      xkp1[i] = xk[i] + ak * pk[i];
      rkp1[i] = rk[i] - ak * Apk[i];
      dkp1 += rkp1[i] * rkp1[i];
    }
    if (dkp1 / numParams < eps || k == kmax)
    {
      converged = true;
      break;
    }
    else
    {
      Real bk = dkp1 / dk;
      for (int i = 0; i < numParams; i++)
      {
        pkp1[i] = rkp1[i] + bk * pk[i];
        pk[i]   = pkp1[i];
        rk[i]   = rkp1[i];
        xk[i]   = xkp1[i];
        dk      = dkp1;
      }
      k++;
    }
  }

  solution = xkp1;
  return k;
}

} // namespace qmcplusplus

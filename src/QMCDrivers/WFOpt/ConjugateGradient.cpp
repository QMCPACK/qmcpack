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
    : threshold_(threshold), regularization_(regularization), conjugate_gradient_timer_(createGlobalTimer("ConjugateGradient::run", timer_level_fine))
{}

int ConjugateGradient::run(QMCCostFunctionBase& optTarget, const std::vector<Real>& bvec, std::vector<Real>& solution)
{
  ScopedTimer local(conjugate_gradient_timer_);
  int numParams = optTarget.getNumParams();
  int kmax      = numParams;
  int k         = 0;

  std::vector<Real> rk(numParams, 0);
  std::vector<Real> rkp1(numParams, 0);
  std::vector<Real> pk(numParams, 0);
  std::vector<Real> pkp1(numParams, 0);
  std::vector<Real> xk(numParams, 0);
  std::vector<Real> xkp1(numParams, 0);
  std::vector<Real> Apk(numParams, 0);

  //Basic implementation of conjugate gradient algorithm...coming directly from wikipedia
  Real dk = 0;
  //initial guess is zero, so A*x_0 = 0
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
    //This function avoids building the full Aij matrix by calculating the matrix vector
    //product of the A matrix * the current solution vector
    //optTarget needs to implement A*p and return in Apk
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

  //store data for potential use in nonlinear rescale
  Axsol_ = std::move(Apk);
  xsol_  = std::move(xkp1);

  //return solution
  solution = xsol_;
  return k;
}

void ConjugateGradient::getNonLinearRange(int& first, int& last, const QMCCostFunctionBase& optTarget) const
{
  std::vector<int> types;
  optTarget.getParameterTypes(types);
  first = 0;
  last  = types.size();
  //assume all non-linear coeffs are together.
  if (types[0] == optimize::LINEAR_P)
  {
    int i(0);
    while (i < types.size())
    {
      if (types[i] == optimize::LINEAR_P)
        first = i;
      i++;
    }
    first++;
  }
  else
  {
    int i(types.size() - 1);
    while (i >= 0)
    {
      if (types[i] == optimize::LINEAR_P)
        last = i;
      i--;
    }
  }
  //     returns the number of non-linear parameters.
  //    app_log()<<"line params: "<<first<<" "<<last<< std::endl;
}

ConjugateGradient::Real ConjugateGradient::getNonLinearRescale(const QMCCostFunctionBase& optTarget) const
{
  int first(0), last(0);
  getNonLinearRange(first, last, optTarget);
  if (first == last)
    return 1.0;
  Real rescale(1.0);
  Real xi(0.5);
  Real D(0.0);
  for (int i = first; i < last; i++)
    D += Axsol_[i] * xsol_[i];
  rescale = (1 - xi) * D / ((1 - xi) + xi * std::sqrt(1 + D));
  rescale = 1.0 / (1.0 - rescale);
  //     app_log()<<"rescale: "<<rescale<< std::endl;
  return rescale;
}

} // namespace qmcplusplus

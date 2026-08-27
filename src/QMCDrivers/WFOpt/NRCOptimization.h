//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NRC_OPTIMIZATION_H
#define QMCPLUSPLUS_NRC_OPTIMIZATION_H

#include <functional>
#include "Configuration.h"

template<class T>
struct NRCOptimization
{
  using Return_t = T;

  /** number of line iteration for Brent method */
  int ITMAX;

  /** maximum CG mixture for line minimization
   *
   * y'=y + Lambda*cg where Lambda < LambdaMax
   */
  Return_t LambdaMax;
  Return_t quadstep;
  Return_t quadoffset;
  Return_t largeQuarticStep;
  Return_t stepsize;
  bool validFuncVal;
  //if tol is absolute, not a percent
  bool AbsFuncTol;

  int current_step;

  NRCOptimization();

  virtual ~NRCOptimization() {}


  Return_t Lambda;
  Return_t ZEPS, CGOLD, TOL, GLIMIT, TINY, GOLD;

  // Returns the number of real roots
  inline int CubicFormula(double a, double b, double c, double d, double& x1, double& x2, double& x3);

  inline Return_t QuarticMinimum(std::vector<Return_t>& coefs);


  bool lineoptimization3(const std::function<Return_t(Return_t)>& evalCost, int points, Return_t& zeroCost);

  bool lineoptimization2(const std::function<Return_t(Return_t)>& evalCost, Return_t maxStep = 1e9);

  inline void shift(Return_t& a, Return_t& b, Return_t& c, Return_t d)
  {
    a = b;
    b = c;
    c = d;
  }

  T brentNRC(const std::function<Return_t(Return_t)>& evalCost, Return_t ax, Return_t bx, Return_t cx, Return_t& xmin);

  bool mnbrakNRC(const std::function<Return_t(Return_t)>& evalCost,
                 Return_t& ax,
                 Return_t& bx,
                 Return_t& cx,
                 Return_t& fa,
                 Return_t& fb,
                 Return_t& fc,
                 Return_t maxStep = 1e9);
};

extern template struct NRCOptimization<qmcplusplus::QMCTraits::RealType>;

#endif

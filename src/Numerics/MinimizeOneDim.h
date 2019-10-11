//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MINIMIZE_ONED_H
#define QMCPLUSPLUS_MINIMIZE_ONED_H

#include <algorithm>
#include <stdexcept>
#include <tuple>

// Storage for bracketed minimum.

template<typename T>
struct Bracket_min_t
{
  T a;
  T b;
  T c;
  bool success;

  Bracket_min_t(T a1, T b1, T c1, bool okay = true) : a(a1), b(b1), c(c1), success(okay) {}
};


// Minimize a function in one dimension
// Bracket a minimum in preparation for minimization

// If 'bound_max' is a positive number, the search range is bounded between zero and 'bound_max'.
// If the search falls outside that range, the function returns with bracket.success set to 'false',
// and the position in bracket.a.   This means the minimum occurs at the edge of the boundary, and
// there is no need to call 'find_minimum' (nor should it be called).


template<class F, typename T>
Bracket_min_t<T> bracket_minimum(const F& f, T initial_value, T bound_max = -1.0)
{
  T xa = initial_value;
  T fa = f(xa);

  T xb = xa + 0.005;
  T fb = f(xb);

  // swap a and b
  if (fb > fa)
  {
    std::swap(xa, xb);
    std::swap(fa, fb);
  }

  bool check_bound = false;
  if (bound_max > 0.0)
  {
    check_bound = true;
  }
  T best_val = xb;

  T delx = 1.61 * (xb - xa);
  T xd   = xb + delx;
  T fd   = f(xd);

  int cnt = 0;
  while (fb > fd)
  {
    T xtmp = xb;
    T ftmp = fb;
    xb     = xd;
    fb     = fd;
    xa     = xtmp;
    fa     = ftmp;
    xd += delx;
    if (check_bound && (xd < 0.0 || xd > bound_max))
    {
      // minimum occurs at the boundary of the range
      return Bracket_min_t<T>(best_val, 0.0, 0.0, false);
    }


    fd = f(xd);

    if (cnt == 50)
    {
      delx *= 5;
    }
    if (cnt == 100)
    {
      delx *= 5;
    }
    if (cnt == 1000)
    {
      delx *= 5;
    }
    if (cnt == 10000)
    {
      delx *= 5;
    }
    cnt++;
    if (cnt == 100000)
    {
      throw std::runtime_error("Failed to bracket minimum");
    }
  }
  if (xa > xd)
    std::swap(xa, xd);
  return Bracket_min_t<T>(xa, xb, xd);
}

// Use a golden-section search

// Returns a pair with the location of the minimum and the value of the function.

template<class F, typename T>
std::pair<T, T> find_minimum(const F& f, Bracket_min_t<T>& bracket)
{
  // assert(bracket.success == true);
  T xa = bracket.a;
  T xb = bracket.b;
  T xd = bracket.c;
  T fb = f(xb);
  T xc = xb + 0.4 * (xd - xb);
  T fc = f(xc);

  T tol = 1e-5;
  while (std::abs(xa - xd) > tol * (std::abs(xb) + std::abs(xc)))
  {
    if (fb > fc)
    {
      xa = xb;
      xb = xa + 0.4 * (xd - xa);
      fb = f(xb);
      xc = xa + 0.6 * (xd - xa);
      fc = f(xc);
    }
    else
    {
      xd = xc;
      xb = xa + 0.4 * (xd - xa);
      fb = f(xb);
      xc = xa + 0.6 * (xd - xa);
      fc = f(xc);
    }
  }
  T final_value;
  T final_x;
  if (fb < fc)
  {
    final_x = xb;
  }
  else
  {
    final_x = xc;
  }
  final_value = f(final_x);
  return std::pair<T, T>(final_x, final_value);
}

#endif

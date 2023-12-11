//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineEval.hpp
 *
 * Master header file to define MultiBspline and MultiBspline1D evaluation functions
 *
 * The evaluation functions in MultiBspline and MultiBspline1D are memory
 * bandwidth (BW) bound. Optimizations towards maximizing BW usage is always
 * needed. For this reason, with SIMD, memory alignment is required in order
 * to saturate the BW. The number of splines must be a multiple of aligned
 * size. The result vectors must be in Structure-of-Array data layout with
 * their starting address correctly aligned.
 *
 */
#ifndef SPLINE2_MULTIEINSPLINE_EVAL_HPP
#define SPLINE2_MULTIEINSPLINE_EVAL_HPP

#include <algorithm>
#include "spline2/bspline_traits.hpp"
#include "spline2/MultiBsplineEval_helper.hpp"

// select evaluation functions based on the architecture.
///include evaluate_v_impl
#include "spline2/MultiBsplineValue.hpp"
#include "spline2/MultiBsplineVGLH.hpp"

///include evaluate_vghgh_impl
#include "spline2/MultiBsplineVGHGH.hpp"

namespace spline2
{
/// evaluate values optionally in the range [first,last)
template<typename SPLINET, typename PT, typename VT>
inline void evaluate3d(const SPLINET& spline, const PT& r, VT& psi)
{
  evaluate_v_impl(spline, r[0], r[1], r[2], psi.data(), 0, psi.size());
}

template<typename SPLINET, typename PT, typename VT>
inline void evaluate3d(const SPLINET& spline, const PT& r, VT& psi, int first, int last)
{
  evaluate_v_impl(spline, r[0], r[1], r[2], psi.data() + first, first, last);
}

/// evaluate values, gradients, laplacians optionally in the range [first,last)
template<typename SPLINET, typename PT, typename VT, typename GT, typename LT>
inline void evaluate3d_vgl(const SPLINET& spline, const PT& r, VT& psi, GT& grad, LT& lap)
{
  evaluate_vgl_impl(spline, r[0], r[1], r[2], psi.data(), grad.data(), lap.data(), psi.size(), 0, psi.size());
}

template<typename SPLINET, typename PT, typename VT, typename GT, typename LT>
inline void evaluate3d_vgl(const SPLINET& spline, const PT& r, VT& psi, GT& grad, LT& lap, int first, int last)
{
  evaluate_vgl_impl(spline, r[0], r[1], r[2], psi.data() + first, grad.data() + first, lap.data() + first, psi.size(),
                    first, last);
}

/// evaluate values, gradients, hessians optionally in the range [first,last)
template<typename SPLINET, typename PT, typename VT, typename GT, typename HT>
inline void evaluate3d_vgh(const SPLINET& spline, const PT& r, VT& psi, GT& grad, HT& hess)
{
  evaluate_vgh_impl(spline, r[0], r[1], r[2], psi.data(), grad.data(), hess.data(), psi.size(), 0, psi.size());
}

template<typename SPLINET, typename PT, typename VT, typename GT, typename HT>
inline void evaluate3d_vgh(const SPLINET& spline, const PT& r, VT& psi, GT& grad, HT& hess, int first, int last)
{
  evaluate_vgh_impl(spline, r[0], r[1], r[2], psi.data() + first, grad.data() + first, hess.data() + first, psi.size(),
                    first, last);
}

template<typename SPLINET, typename PT, typename VT, typename GT, typename HT, typename GHT>
inline void evaluate3d_vghgh(const SPLINET& spline, const PT& r, VT& psi, GT& grad, HT& hess, GHT& ghess)
{
  evaluate_vghgh_impl(spline, r[0], r[1], r[2], psi.data(), grad.data(), hess.data(), ghess.data(), psi.size(), 0,
                      psi.size());
}

template<typename SPLINET, typename PT, typename VT, typename GT, typename HT, typename GHT>
inline void evaluate3d_vghgh(const SPLINET& spline,
                                    const PT& r,
                                    VT& psi,
                                    GT& grad,
                                    HT& hess,
                                    GHT& ghess,
                                    int first,
                                    int last)
{
  evaluate_vghgh_impl(spline, r[0], r[1], r[2], psi.data() + first, grad.data() + first, hess.data() + first,
                      ghess.data() + first, psi.size(), first, last);
}

} // namespace spline2
#endif

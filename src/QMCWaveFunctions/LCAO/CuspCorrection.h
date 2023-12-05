//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file CuspCorrection.h
  * @brief Corrections to electron-nucleus cusp for all-electron molecular calculations.
  */

#ifndef QMCPLUSPLUS_CUSPCORRECTION_H
#define QMCPLUSPLUS_CUSPCORRECTION_H

#include <cmath>
#include "Configuration.h"

namespace qmcplusplus
{
/**
  * @brief Cusp correction parameters
  *
  * From "Scheme for adding electron-nuclear cusps to Gaussian orbitals"  Ma, Towler, Drummond, and Needs
  *  JCP 122, 224322 (2005)
  *
  * Equations 7 and 8 in the paper define the correction.  These are the parameters in those equations.
  */

struct CuspCorrectionParameters
{
  using ValueType = QMCTraits::ValueType;
  using RealType  = QMCTraits::RealType;

  /// The cutoff radius
  RealType Rc;

  /// A shift to keep correction to a single sign
  RealType C;

  /// The sign of the wavefunction at the nucleus
  RealType sg;

  /// The coefficients of the polynomial \f$p(r)\f$ in Eq 8
  TinyVector<ValueType, 5> alpha;

  /// Flag to indicate the correction should be recalculated
  int redo;

  CuspCorrectionParameters() : Rc(0.0), C(0.0), sg(1.0), alpha(0.0), redo(0) {}
};

/// Formulas for applying the cusp correction

class CuspCorrection
{
  using RealType = QMCTraits::RealType;

public:
  inline RealType Rr(RealType r) const { return cparam.sg * std::exp(pr(r)); }

  inline RealType pr(RealType r) const
  {
    auto& alpha = cparam.alpha;
    return alpha[0] + alpha[1] * r + alpha[2] * r * r + alpha[3] * r * r * r + alpha[4] * r * r * r * r;
  }

  inline RealType dpr(RealType r) const
  {
    auto& alpha = cparam.alpha;
    return alpha[1] + 2.0 * alpha[2] * r + 3.0 * alpha[3] * r * r + 4.0 * alpha[4] * r * r * r;
  }

  inline RealType d2pr(RealType r) const
  {
    auto& alpha = cparam.alpha;
    return 2.0 * alpha[2] + 6.0 * alpha[3] * r + 12.0 * alpha[4] * r * r;
  }

  CuspCorrection(const CuspCorrectionParameters& param) : cparam(param) {}

  CuspCorrectionParameters cparam;
};
} // namespace qmcplusplus

#endif

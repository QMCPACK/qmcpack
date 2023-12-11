//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NANGUARD_H
#define QMCPLUSPLUS_NANGUARD_H

#include "Configuration.h"

namespace qmcplusplus
{

class NaNguard
{
  using PsiValue = QMCTraits::QTFull::ValueType;
  using GradType = QMCTraits::GradType;

public:
  /** check if ratio is NaN and throw an error if yes.
   * @param ratio psi ratio to be checked
   * @param message printout to indicate what the issue is.
   */
  static void checkOneParticleRatio(const PsiValue& ratio, const std::string_view info);

  /** check if any gradient component (x,y,z) is NaN and throw an error if yes.
   * @param grads gradients to be checked
   * @param message printout to indicate what the issue is.
   */
  static void checkOneParticleGradients(const GradType& grads, const std::string_view info);
};
} // namespace qmcplusplus
#endif

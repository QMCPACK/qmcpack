//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_WAVEFUNCTIONTESTERBATCHEDINPUT_H
#define QMCPLUSPLUS_WAVEFUNCTIONTESTERBATCHEDINPUT_H

#include "Configuration.h"

namespace qmcplusplus
{
class WaveFunctionTesterBatchedInput
{
public:
  using RealType = QMCTraits::RealType;

  /** Parse supported batched tester controls and reject legacy-only diagnostics. */
  void readXML(xmlNodePtr xml_input);

  /** Return whether the finite-difference gradient and Laplacian diagnostic is enabled. */
  bool get_run_basic() const { return run_basic_; }
  /** Return whether the ratio consistency diagnostic is enabled. */
  bool get_run_ratio() const { return run_ratio_; }
  /** Return the finite-difference displacement, or zero when the default should be used. */
  RealType get_delta() const { return delta_; }
  /** Return the relative error tolerance, or zero when the default should be used. */
  RealType get_tolerance() const { return tolerance_; }

private:
  bool run_basic_ = true;
  bool run_ratio_ = false;
  RealType delta_ = 0.0;
  RealType tolerance_ = 0.0;
};
} // namespace qmcplusplus

#endif
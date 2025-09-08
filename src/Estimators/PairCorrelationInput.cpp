//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/PairCorrEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "PairCorrelationInput.h"
#include "EstimatorInput.h"

namespace qmcplusplus
{

PairCorrelationInput::PairCorrelationInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(name_, "name");
  setIfInInput(type_, "type");
  explicit_set_nbins_ = setIfInInput(nbins_, "num_bin");
  explicit_set_rmax_  = setIfInInput(rmax_, "rmax");
  explicit_set_delta_ = setIfInInput(delta_, "dr");
  setIfInInput(debug_, "debug");
  setIfInInput(sources_, "sources");
  setIfInInput(target_, "target");
}

} // namespace qmcplusplus

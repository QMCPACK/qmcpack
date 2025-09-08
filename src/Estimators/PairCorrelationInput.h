//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/PairCorrEstimator.h
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_PAIRCORRELATIONINPUT_H
#define QMCPLUSPLUS_PAIRCORRELATIONINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{

class PairCorrelationEstimator;

class PairCorrelationInput
{
public:
  using Consumer = PairCorrelationEstimator;
  using Real     = QMCTraits::FullPrecRealType;

  static constexpr std::string_view type_tag{"PairCorrelation"};
  class PairCorrelationInputSection : public InputSection
  {
  public:
    PairCorrelationInputSection()
    {
      section_name            = type_tag;
      section_name_alternates = {"gofr"};
      attributes              = {"type", "name", "num_bin", "rmax", "dr", "debug", "sources"};
      strings                 = {"type", "name"};
      multi_strings           = {"sources"};
      reals                   = {"dr", "rmax"};
      integers                = {"num_bin"};
      bools                   = {"debug"};
    }
  };

  PairCorrelationInput(xmlNodePtr);

private:
  PairCorrelationInputSection input_section_;

  std::string name_{type_tag};
  std::string type_{type_tag};
  std::vector<std::string> sources_{"e"};
  Real rmax_{10};
  bool explicit_set_rmax_{false};
  Real delta_{0.5};
  bool explicit_set_delta_{false};
  int nbins_{20};
  bool explicit_set_nbins_{false};
  bool debug_{false};

public:
  std::string get_name() const { return name_; }
  std::string get_type() const { return type_; }
  const std::vector<std::string>& get_sources() const { return sources_; }
  Real get_rmax() const { return rmax_; }
  Real get_delta() const { return delta_; }
  int get_nbins() const { return nbins_; }
  bool get_debug() const { return debug_; }
  bool get_explicit_set_rmax() const { return explicit_set_rmax_; }
  bool get_explicit_set_delta() const { return explicit_set_delta_; };
  bool get_explicit_set_nbins() const { return explicit_set_nbins_; };
};

} // namespace qmcplusplus
#endif

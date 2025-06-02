//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_SHOVERLAPINPUT_H
#define QMCPLUSPLUS_SHOVERLAPINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{

class SelfHealingOverlap;

/** Native representation for Self-Healing Overlap Estimator inputs
 */
class SelfHealingOverlapInput
{
public:
  static constexpr std::string_view type_tag{"SelfHealingOverlap"};

  using Consumer = SelfHealingOverlap;
  using Real     = QMCTraits::RealType;

  class SelfHealingOverlapInputSection : public InputSection
  {
  public:
    // clang-format: off
    SelfHealingOverlapInputSection()
    {
      section_name   = type_tag;
      attributes     = {"type", "name", "param_deriv"};
      strings        = {"type", "name"};
      bools          = {"param_deriv"};
      default_values = {{"param_deriv", false}};
    }
    // clang-format: on
  };

  std::string get_name() const { return name_; }
  std::string get_type() const { return type_; }

  SelfHealingOverlapInput(xmlNodePtr cur);
  /** default copy constructor
   *  This is required due to MDI being part of a variant used as a vector element.
   */
  SelfHealingOverlapInput(const SelfHealingOverlapInput&) = default;

  // SelfHealingOverlap violates the encapsulation of SelfHealingOverlapInput
  SelfHealingOverlapInputSection input_section_;

private:
  std::string name_{type_tag};
  std::string type_{type_tag};
};

} // namespace qmcplusplus
#endif /* SHOVERLAPINPUT_H */

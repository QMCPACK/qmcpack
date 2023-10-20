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
  using Consumer = SelfHealingOverlap;
  using Real = QMCTraits::RealType;

  class SelfHealingOverlapInputSection : public InputSection
  {
  public:
    // clang-format: off
    SelfHealingOverlapInputSection()
    {
      section_name   = "SelfHealingOverlap";
      attributes     = {"type", "name"};
      strings        = {"type", "name"};
      default_values = {{"type", std::string("sh_overlap")},{"name", std::string("sh_overlap")}};
    }
    // clang-format: on
  };

  SelfHealingOverlapInput(xmlNodePtr cur) { input_section_.readXML(cur); }
  /** default copy constructor
   *  This is required due to MDI being part of a variant used as a vector element.
   */
  SelfHealingOverlapInput(const SelfHealingOverlapInput&) = default;
  SelfHealingOverlapInputSection input_section_;
};

} // namespace qmcplusplus
#endif /* SHOVERLAPINPUT_H */

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_DIPOLEMOMENTINPUT_H
#define QMCPLUSPLUS_DIPOLEMOMENTINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{

class DipoleMoment;

/** Native representation for Self-Healing Overlap Estimator inputs
 */
class DipoleMomentInput
{
public:
  using Consumer = DipoleMoment;
  using Real = QMCTraits::RealType;

  class DipoleMomentInputSection : public InputSection
  {
  public:
    // clang-format: off
    DipoleMomentInputSection()
    {
      section_name   = "DipoleMoment";
      attributes     = {"type", "name", "ions"};
      strings        = {"type", "name", "ions"};
      default_values = {{"type"     , std::string("dipolemoment")},
                        {"name"     , std::string("dipolemoment")},
                        {"ions"     , std::string("ion0"        )}};
    }
    // clang-format: on
  };

  DipoleMomentInput(xmlNodePtr cur) { input_section_.readXML(cur); }
  /** default copy constructor
   *  This is required due to MDI being part of a variant used as a vector element.
   */
  DipoleMomentInput(const DipoleMomentInput&) = default;
  DipoleMomentInputSection input_section_;
};

} // namespace qmcplusplus
#endif /* DIPOLEMOMENTINPUT_H */

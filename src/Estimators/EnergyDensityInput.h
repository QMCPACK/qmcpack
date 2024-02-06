//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/EnergyDensityEstimator.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ENERGY_DENSITY_INPUT_H
#define QMCPLUSPLUS_ENERGY_DENSITY_INPUT_H

#include "InputSection.h"
#include "ReferencePointsInput.h"
#include "SpaceGridInput.h"
#include "type_traits/template_types.hpp"

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class EnergyDensityTests;
}

class NEEnergyDensityEstimator;
template<typename REAL>
class NESpaceGrid;

/** EnergyDensity has two other XML input reading objects that it delegates to.
 *  I don't think these need to be handled with a variant vector because there
 *  is no expectation it be expanded and they in general occur together and not as options.
 */
// using EnergyDensityDelegate =
//     std::variant<std::monostate, ReferencePointsInput, SpaceGridInput >;
// using EnergyDensityDelegates = std::vector<EnergyDensityDelegate>;

/** Native representation for DensityMatrices1B Estimator's inputs
 */
class EnergyDensityInput
{
public:
  using Consumer = NEEnergyDensityEstimator;

  class EnergyDensityInputSection : public InputSection
  {
  public:
    EnergyDensityInputSection()
    {
      section_name = "EnergyDensity";
      attributes   = {"name", "dynamic", "static", "ion_points", "type"};
      parameters   = {"reference_points", "spacegrid"};
      strings      = {"name", "type", "dynamic", "static"};
      bools        = {"ion_points"};
      delegates    = {"reference_points", "spacegrid"};
      multiple     = {"spacegrid"};
      required     = {"spacegrid"};
      registerDelegate("reference_points", makeReferencePointsInput);
      registerDelegate("spacegrid", makeSpaceGridInput);
    }
    /** Here the delegate input object is registered */
    EnergyDensityInputSection(const EnergyDensityInputSection&) = default;
  };

  using Real = QMCTraits::RealType;

  EnergyDensityInput(const EnergyDensityInput&) = default;
  EnergyDensityInput(xmlNodePtr cur);

  const std::string& get_name() const { return name_; }
  const std::string& get_dynamic() const { return dynamic_; }
  const std::string& get_static() const { return static_; }
  ReferencePointsInput get_ref_points_input() const { return ref_points_input_; }
  std::vector<SpaceGridInput> get_space_grid_inputs() const;
  const bool& get_ion_points() const { return ion_points_; }
private:
  std::string name_;
  std::string dynamic_;
  std::string static_;
  bool ion_points_{false};
  EnergyDensityInputSection input_section_;
  ReferencePointsInput ref_points_input_;
};

} // namespace qmcplusplus

#endif

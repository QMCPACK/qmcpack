//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: SpaceGrid.h & SpaceGrid.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPACEGRID_INPUT_H
#define QMCPLUSPLUS_SPACEGRID_INPUT_H

#include <set>
#include <unordered_map>

#include <Configuration.h>
#include "InputSection.h"
#include "EstimatorInput.h"
#include "ParseGridInput.hpp" //for AxisGrid

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class EnergyDensityTests;
}

class SpaceGrid;

class SpaceGridInput
{
public:
  using Consumer = SpaceGrid;
  using Real     = double;

  enum class CoordForm
  {
    CARTESIAN = 0,
    CYLINDRICAL,
    SPHERICAL
  };

  inline static const std::unordered_map<std::string, std::any> lookup_input_enum_value{{"coord-cartesian",
                                                                                         CoordForm::CARTESIAN},
                                                                                        {"coord-cylindrical",
                                                                                         CoordForm::CYLINDRICAL},
                                                                                        {"coord-spherical",
                                                                                         CoordForm::SPHERICAL}};

  using LabelSet = std::vector<std::string_view>;
  // legal labels for each coordinate type.  These are effectively enums
  inline static const LabelSet ax_cartesian{"x", "y", "z"};
  inline static const LabelSet ax_cylindrical{"r", "phi", "z"};
  inline static const LabelSet ax_spherical{"r", "phi", "theta"};
  inline static const std::unordered_map<CoordForm, LabelSet> axes_label_sets{{CoordForm::CARTESIAN, ax_cartesian},
                                                                               {CoordForm::CYLINDRICAL, ax_cylindrical},
                                                                               {CoordForm::SPHERICAL, ax_spherical}};

  class SpaceGridAxisInput
  {
    class SpaceGridAxisInputSection : public InputSection
    {
    public:
      SpaceGridAxisInputSection()
      {
        section_name = "axis";
        attributes   = {"label", "grid", "p1", "p2", "scale"};
        strings      = {"label", "p1", "p2"};
        custom = {"grid"}, reals = {"scale"};
        required = {"label", "p1"};
      }
      void setFromStreamCustom(const std::string& ename, const std::string& name, std::istringstream& svalue) override;
      SpaceGridAxisInputSection(const SpaceGridAxisInputSection& sgais) = default;
    };

  public:
    SpaceGridAxisInput(xmlNodePtr cur);
    SpaceGridAxisInput(const SpaceGridAxisInput& sgai) = default;

    static std::any makeAxis(xmlNodePtr cur, std::string& value_label)
    {
      SpaceGridAxisInput space_grid_axis{cur};

      value_label = "axis";
      return space_grid_axis;
    }

    const SpaceGridAxisInputSection& get_input() { return input_section_; }

    std::string get_label() const { return label_; }
    Real get_scale() const { return scale_; }
    std::string get_p1() const { return p1_; }
    std::string get_p2() const { return p2_; }
    AxisGrid<Real> get_grid() const { return grid_; }

  private:
    SpaceGridAxisInputSection input_section_;
    std::string label_ = "";
    Real scale_        = 1.0;
    std::string p1_    = "";
    std::string p2_    = "zero";
    AxisGrid<Real> grid_;
  };

  class SpaceGridOriginInput
  {
    class SpaceGridOriginInputSection : public InputSection
    {
    public:
      SpaceGridOriginInputSection()
      {
        section_name = "origin";
        attributes   = {"p1", "p2", "fraction"};
        required     = {"p1"};
        strings      = {"p1", "p2"};
        reals        = {"fraction"};
      }
      SpaceGridOriginInputSection(const SpaceGridOriginInputSection& sgois) = default;
    };

  public:
    SpaceGridOriginInput(xmlNodePtr cur);

    static std::any makeOrigin(xmlNodePtr cur, std::string& value_label)
    {
      SpaceGridOriginInput space_grid_origin{cur};
      value_label = "origin";
      return space_grid_origin;
    }
    const std::string& get_p1() const { return p1_; }
    const std::string& get_p2() const { return p2_; }
    const Real get_fraction() const { return fraction_; }
  private:
    SpaceGridOriginInputSection input_section_;
    std::string p1_{"zero"};
    std::string p2_{""};
    Real fraction_{0.0};
  };

  class SpaceGridInputSection : public InputSection
  {
  public:
    SpaceGridInputSection()
    {
      section_name = "SpaceGrid";
      attributes   = {"coord"};
      enums        = {"coord"};
      delegates    = {"origin", "axis"};
      multiple     = {"axis"};
      registerDelegate("origin", SpaceGridOriginInput::makeOrigin);
      registerDelegate("axis", SpaceGridAxisInput::makeAxis);
    }
    std::any assignAnyEnum(const std::string& name) const override;
    void checkParticularValidity() override;
    SpaceGridInputSection(const SpaceGridInputSection& sgis) = default;
  };

  SpaceGridInput(xmlNodePtr cur);
  SpaceGridInput(const SpaceGridInput& sgi) = default;

  CoordForm get_coord_form() const { return coord_form_; }
  bool isPeriodic() const { return !(coord_form_ == CoordForm::CYLINDRICAL || coord_form_ == CoordForm::SPHERICAL); }
  const std::array<std::string, OHMMS_DIM>& get_axis_p1s() const { return axis_p1s_; }
  const std::array<std::string, OHMMS_DIM>& get_axis_p2s() const { return axis_p2s_; }

  const std::array<Real, OHMMS_DIM>& get_axis_scales() const { return axis_scales_; }
  const std::array<std::string, OHMMS_DIM>& get_axis_labels() const { return axis_labels_; }
  const std::array<AxisGrid<Real>, OHMMS_DIM>& get_axis_grids() const { return axis_grids_; }
  const std::string& get_origin_p1() const { return origin_p1_; }
  const std::string& get_origin_p2() const { return origin_p2_; }
  Real get_origin_fraction() const { return origin_fraction_; }
  /** axes_label_set accessor, avoids a bunch of switch statements
   *  at must be used because std::unordered_map::operator[] can't return a const reference
   */
  const LabelSet& get_axes_label_set() const { return axes_label_sets.at(coord_form_); }

private:
  void checkAxes(std::vector<std::any>& axes);
  void checkGrids();

  SpaceGridInputSection input_section_;
  CoordForm coord_form_;
  // origin in optional so this is required.
  std::string origin_p1_{"zero"};
  std::string origin_p2_{""};
  Real origin_fraction_{0.0};
  std::array<std::string, OHMMS_DIM> axis_labels_;
  std::array<std::string, OHMMS_DIM> axis_p1s_;
  std::array<std::string, OHMMS_DIM> axis_p2s_;
  std::array<Real, OHMMS_DIM> axis_scales_;
  std::array<AxisGrid<Real>, OHMMS_DIM> axis_grids_;
};

/** factory function for a SpaceGridInput
 *  \param[in]  input node for SpaceGridInput
 *  \param[out] value label returned to caller
 */
std::any makeSpaceGridInput(xmlNodePtr, std::string& value_label);

} // namespace qmcplusplus
#endif
